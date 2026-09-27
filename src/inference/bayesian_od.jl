# Turing.jl Probabilistic Model for Bayesian Orbit Determination & Parameter Estimation

"""
    bayesian_orbit_determination_model(t_obs, range_obs, doppler_obs, gs, prior_u0_mean, prior_u0_std, opts)

Turing `@model` macro defining Bayesian posterior over initial state u0 = [x, y, z, vx, vy, vz] and measurement noise parameters.
"""
@model function bayesian_orbit_determination_model(t_obs, range_obs, doppler_obs, gs::GroundStation, prior_u0_mean::Vector{Float64}, prior_u0_std::Vector{Float64}, opts::OrbitPropagatorOptions)
    # Sample initial state parameters with concrete type inference for Turing AD
    u1 ~ Normal(prior_u0_mean[1], prior_u0_std[1])
    u2 ~ Normal(prior_u0_mean[2], prior_u0_std[2])
    u3 ~ Normal(prior_u0_mean[3], prior_u0_std[3])
    u4 ~ Normal(prior_u0_mean[4], prior_u0_std[4])
    u5 ~ Normal(prior_u0_mean[5], prior_u0_std[5])
    u6 ~ Normal(prior_u0_mean[6], prior_u0_std[6])

    T_elem = typeof(u1)
    u0 = T_elem[u1, u2, u3, u4, u5, u6]

    # Sensor measurement noise specifications (10m range, 1mm/s doppler)
    sigma_range = 0.010   # km
    sigma_doppler = 1e-4  # km/s

    # Propagate orbit dynamics ODE
    tspan = (t_obs[1], t_obs[end])
    sol = propagate_orbit(u0, tspan, t_obs, opts)

    # Likelihood integration over observation epochs
    for (i, t) in enumerate(t_obs)
        r_sat = SVector{3}(sol.u[i][1], sol.u[i][2], sol.u[i][3])
        v_sat = SVector{3}(sol.u[i][4], sol.u[i][5], sol.u[i][6])
        meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t)

        range_obs[i] ~ Normal(meas.range, sigma_range)
        doppler_obs[i] ~ Normal(meas.range_rate, sigma_doppler)
    end
end

"""
    fit_bayesian_od(t_obs, range_obs, doppler_obs, gs, prior_mean, prior_std, opts; n_samples=500, n_adapt=250, target_accept=0.90)

Run NUTS MCMC sampling for orbit determination state and noise posterior inference.
Uses target_accept=0.90 for smooth ODE Hamiltonian dynamics.
"""
function fit_bayesian_od(t_obs::Vector{Float64}, range_obs::Vector{Float64}, doppler_obs::Vector{Float64}, gs::GroundStation, prior_mean::Vector{Float64}, prior_std::Vector{Float64}, opts::OrbitPropagatorOptions; n_samples::Int=500, n_adapt::Int=250, target_accept::Float64=0.90)
    model = bayesian_orbit_determination_model(t_obs, range_obs, doppler_obs, gs, prior_mean, prior_std, opts)
    chain = sample(model, NUTS(n_adapt, target_accept), n_samples, progress=false)
    return chain
end
