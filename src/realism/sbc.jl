# Simulation-Based Calibration (SBC) Engine for Orbit Determination

"""
    compute_sbc_rank_statistics(n_simulations, gs, prior_mean, prior_std, opts; n_samples=300)

Run Simulation-Based Calibration (SBC) over `n_simulations` trials:
1. Draw ground truth parameter u_true ~ Prior
2. Generate synthetic tracking data y_sim ~ Likelihood(u_true)
3. Fit Bayesian posterior using MCMC
4. Compute rank of u_true within MCMC posterior draws

If Bayesian OD algorithm is well-calibrated, rank statistics follow a uniform distribution.
"""
function compute_sbc_rank_statistics(n_simulations::Int, gs::GroundStation, prior_mean::Vector{Float64}, prior_std::Vector{Float64}, opts::OrbitPropagatorOptions; n_samples::Int=300, n_adapt::Int=100)
    dim = length(prior_mean)
    ranks = Matrix{Int}(undef, n_simulations, dim)

    t_obs = collect(range(0.0, stop=3600.0, length=20))
    sigma_r_true = 0.01 # 10m range noise
    sigma_d_true = 1e-4 # 1mm/s doppler noise

    for sim in 1:n_simulations
        # 1. Sample ground truth
        u_true = [rand(Normal(prior_mean[j], prior_std[j])) for j in 1:dim]

        # 2. Simulate tracking observations
        sim_data = simulate_tracking_data(u_true, t_obs, gs, sigma_r_true, sigma_d_true; opts=opts)

        # 3. Fit MCMC
        chain = fit_bayesian_od(t_obs, sim_data.range, sim_data.doppler, gs, prior_mean, prior_std, opts; n_samples=n_samples, n_adapt=n_adapt)

        # 4. Extract samples and compute rank
        u_samples = hcat([vec(chain[k]) for k in [:u1, :u2, :u3, :u4, :u5, :u6]]...)
        for j in 1:dim
            ranks[sim, j] = count(x -> x < u_true[j], u_samples[:, j])
        end
    end

    return ranks
end
