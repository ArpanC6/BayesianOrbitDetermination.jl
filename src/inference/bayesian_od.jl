# Turing.jl Probabilistic Model for Bayesian Orbit Determination

"""
    bayesian_orbit_determination_model(t_obs, range_obs, doppler_obs, gs,
                                        prior_u0_mean, prior_u0_std, opts)

Turing `@model` for Bayesian orbit determination.

The posterior is over the initial state
`u0 = [x, y, z, vx, vy, vz]` with Gaussian priors centered at `prior_u0_mean`
and standard deviations `prior_u0_std`. Observation noise is fixed to the
values used in the synthetic data generator (10 m range, 1 mm/s Doppler).
"""
@model function bayesian_orbit_determination_model(
    t_obs, range_obs, doppler_obs, gs::GroundStation,
    prior_u0_mean::Vector{Float64}, prior_u0_std::Vector{Float64},
    opts::OrbitPropagatorOptions,
)
    # Non-centered parameterization (CRITICAL for NUTS to work properly)
    z1 ~ Normal(0.0, 1.0)
    z2 ~ Normal(0.0, 1.0)
    z3 ~ Normal(0.0, 1.0)
    z4 ~ Normal(0.0, 1.0)
    z5 ~ Normal(0.0, 1.0)
    z6 ~ Normal(0.0, 1.0)

    # Transform back to physical space (km and km/s)
    u1 = prior_u0_mean[1] + prior_u0_std[1] * z1
    u2 = prior_u0_mean[2] + prior_u0_std[2] * z2
    u3 = prior_u0_mean[3] + prior_u0_std[3] * z3
    u4 = prior_u0_mean[4] + prior_u0_std[4] * z4
    u5 = prior_u0_mean[5] + prior_u0_std[5] * z5
    u6 = prior_u0_mean[6] + prior_u0_std[6] * z6

    u0 = [u1, u2, u3, u4, u5, u6]

    sigma_range = 0.010   # km
    sigma_doppler = 1e-4  # km/s

    tspan = (t_obs[1], t_obs[end])
    sol = propagate_orbit(u0, tspan, t_obs, opts)

    for (i, t) in enumerate(t_obs)
        r_sat = SVector{3}(sol.u[i][1], sol.u[i][2], sol.u[i][3])
        v_sat = SVector{3}(sol.u[i][4], sol.u[i][5], sol.u[i][6])
        meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t)

        range_obs[i] ~ Normal(meas.range, sigma_range)
        doppler_obs[i] ~ Normal(meas.range_rate, sigma_doppler)
    end
end

"""
    fit_bayesian_od(t_obs, range_obs, doppler_obs, gs, prior_mean, prior_std, opts;
                    n_samples = 500, n_adapt = 250, target_accept = 0.90)

Run NUTS MCMC for the orbit determination posterior.

# Arguments
- `t_obs`, `range_obs`, `doppler_obs`: observation times and measurements
- `gs`         : ground station
- `prior_mean` : prior mean for `u0` (6-vector)
- `prior_std`  : prior standard deviations for `u0` (6-vector)
- `opts`       : propagation options
- `n_samples`  : number of posterior samples
- `n_adapt`    : number of adaptation steps
- `target_accept`: NUTS target acceptance rate

# Returns
An MCMC chain with samples for `z1, ..., z6`.
"""
function fit_bayesian_od(
    t_obs::Vector{Float64}, range_obs::Vector{Float64}, doppler_obs::Vector{Float64},
    gs::GroundStation, prior_mean::Vector{Float64}, prior_std::Vector{Float64},
    opts::OrbitPropagatorOptions;
    n_samples::Int = 500, n_adapt::Int = 250, target_accept::Float64 = 0.90,
)
    model = bayesian_orbit_determination_model(
        t_obs, range_obs, doppler_obs, gs, prior_mean, prior_std, opts,
    )
    # Fixed NUTS syntax (n_adapt is positional, not keyword)
    chain = sample(
        model,
        NUTS(n_adapt, target_accept),
        n_samples;
        progress = false,
    )
    return chain
end