# Tutorial 03: Frequentist (Batch Least Squares / Orekit baseline) vs Bayesian Orbit Determination

using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf

println("--- BayesianOrbitDetermination.jl Tutorial 03 ---")
println("Comparison: Frequentist Nonlinear Least Squares vs Bayesian MCMC Posterior\n")

kepler_true = KeplerianState(7000.0, 0.005, deg2rad(28.5), deg2rad(10.0), deg2rad(15.0), deg2rad(0.0))
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

gs = GroundStation("Kourou", 5.251, -52.805, 0.015)
t_obs = collect(range(0.02, stop=3600.0, length=20))

sigma_r = 0.020 # 20m range noise
sigma_d = 2e-4  # 0.2 m/s doppler noise

opts = OrbitPropagatorOptions(include_j2=true)
data = simulate_tracking_data(u0_true, t_obs, gs, sigma_r, sigma_d; opts=opts)

guess_u0 = u0_true + [1.0, -1.0, 0.5, 0.002, -0.002, 0.001]

println("1. Running Frequentist Batch Least Squares (MLE / Gauss-Newton / Orekit approach)...")
bls_res = fit_batch_least_squares(t_obs, data.range, data.doppler, gs, guess_u0, sigma_r, sigma_d, opts)
u_mle = bls_res.u_estimated
cov_mle = bls_res.covariance

d2_mle = compute_mahalanobis_distance(u0_true, u_mle, cov_mle)

println("\n2. Running Bayesian NUTS MCMC...")
prior_mean = bls_res.u_estimated
prior_std  = [0.5, 0.5, 0.5, 0.001, 0.001, 0.001]
chain = fit_bayesian_od(t_obs, data.range, data.doppler, gs, prior_mean, prior_std, opts;
                         n_samples=1000, n_adapt=500, target_accept=0.80)

z_samples = hcat([vec(chain[k]) for k in [:z1, :z2, :z3, :z4, :z5, :z6]]...)
u_samples = z_samples .* prior_std' .+ prior_mean'

u_bayes_mean = vec(mean(u_samples, dims=1))
cov_bayes = cov(u_samples)

eval_bayes = EvaluateCovarianceRealism(u0_true, u_samples, cov_bayes)

println("\n---------------- COMPARISON RESULTS ----------------")
@printf("Frequentist MLE Position Error: %.4f km | D^2: %.2f\n", norm(u0_true[1:3] - u_mle[1:3]), d2_mle)
@printf("Bayesian Posterior Position Error: %.4f km | D^2: %.2f\n", norm(u0_true[1:3] - u_bayes_mean[1:3]), eval_bayes.mahalanobis_d2)
@printf("Bayesian Chi2 GOF p-value: %.4f (Realistic Covariance: %s)\n", eval_bayes.p_value_gof, string(!eval_bayes.is_overconfident))
println("----------------------------------------------------")