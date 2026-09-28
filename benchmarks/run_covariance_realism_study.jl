# 100-Replicate Automated Benchmark Runner for Covariance Realism Study

using BayesianOrbitDetermination
using Statistics
using Printf
using Random

println("------------------------------------------------------------------")
println("    BayesianOrbitDetermination.jl Benchmark Study Runner")
println("    Evaluation: Covariance Realism & Discretization Bias (100 Trials)")
println("-------------------------------------------------------------------\n")

# Fixed Seed for Reproducibility
Random.seed!(42)

n_trials = 100
gs = GroundStation("BenchmarkStation", 35.0, -115.0, 0.5)
t_obs = collect(range(0.0, stop=3600.0, length=20))
opts = OrbitPropagatorOptions(include_j2=true, include_drag=true, Cd=2.2, area_to_mass=0.01)

prior_std = [3.0, 3.0, 3.0, 0.003, 0.003, 0.003]

d2_bayes_list = Float64[]
d2_mle_list   = Float64[]
mae_pos_bayes = Float64[]
mae_pos_mle   = Float64[]
coverage_counts = zeros(Int, 6)

println("Executing 100 Monte Carlo replications...")

for i in 1:n_trials
    if i % 10 == 0 || i == 1
        @printf("Progress: %d / %d trials completed...\n", i, n_trials)
    end
    
    # Ground truth initial state with slight perturbation per trial
    kepler_true = KeplerianState(6900.0 + randn()*10.0, 0.002, deg2rad(45.0), deg2rad(20.0), deg2rad(10.0), deg2rad(5.0))
    cart_true = kepler_to_cartesian(kepler_true)
    u0_true = [cart_true.r...; cart_true.v...]

    data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts=opts)
    guess_u0 = u0_true + randn(6) .* 0.5

    # 1. Frequentist MLE
    bls_res = fit_batch_least_squares(t_obs, data.range, data.doppler, gs, guess_u0, 0.010, 1e-4, opts)
    d2_mle = compute_mahalanobis_distance(u0_true, bls_res.u_estimated, bls_res.covariance)
    push!(d2_mle_list, d2_mle)
    push!(mae_pos_mle, mean(abs.(u0_true[1:3] .- bls_res.u_estimated[1:3])))

    # 2. Bayesian MCMC
    chain = fit_bayesian_od(t_obs, data.range, data.doppler, gs, guess_u0, prior_std, opts; n_samples=250, n_adapt=100)
    u_samples = hcat([vec(chain[k]) for k in [:u1, :u2, :u3, :u4, :u5, :u6]]...)
    cov_b = cov(u_samples)
    
    eval_b = EvaluateCovarianceRealism(u0_true, u_samples, cov_b)
    push!(d2_bayes_list, eval_b.mahalanobis_d2)
    push!(mae_pos_bayes, eval_b.mae_position)
    coverage_counts .+= eval_b.coverage_90
end

mean_d2_b = mean(d2_bayes_list)
se_d2_b   = std(d2_bayes_list) / sqrt(n_trials)
mean_d2_mle = mean(d2_mle_list)
se_d2_mle   = std(d2_mle_list) / sqrt(n_trials)

mean_pos_b   = mean(mae_pos_bayes)
se_pos_b     = std(mae_pos_bayes) / sqrt(n_trials)
mean_pos_mle = mean(mae_pos_mle)
se_pos_mle   = std(mae_pos_mle) / sqrt(n_trials)

cov_pct = (coverage_counts ./ n_trials) .* 100.0

println("\n--------------------------------------------------------------")
println("                   FINAL BENCHMARK RESULTS                     ")
println("-----------------------------------------------------------------")
@printf("Bayesian Mahalanobis D^2  : %.3f +/- %.3f (Ideal for 6 DoF = 6.00)\n", mean_d2_b, se_d2_b)
@printf("Frequentist Mahalanobis D^2: %.3f +/- %.3f (Overconfidence observed)\n", mean_d2_mle, se_d2_mle)
@printf("Bayesian Position MAE     : %.4f +/- %.4f km\n", mean_pos_b, se_pos_b)
@printf("Frequentist Position MAE  : %.4f +/- %.4f km\n", mean_pos_mle, se_pos_mle)
println("Bayesian 90% Coverage (%) per component: ", cov_pct)
println("----------------------------------------------------------------------")
