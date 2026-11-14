# Turing.jl Probabilistic Model for Bayesian Orbit Determination

"""
    bayesian_orbit_determination_model(t_obs, range_obs, doppler_obs, gs,
                                        prior_u0_mean, prior_u0_std, opts)

Turing `@model` for Bayesian orbit determination.

The posterior is over the initial state
`u0 = [x, y, z, vx, vy, vz]` **at epoch t = 0**, with Gaussian priors centered
at `prior_u0_mean` and standard deviations `prior_u0_std`. Observation times
are seconds from that epoch. Observation noise is Gaussian with standard
deviations `sigma_range` (km) and `sigma_doppler` (km/s); the defaults are
10 m and 1 mm/s.

If `whitening_L` is supplied it must be a lower-triangular factor of a
whitening covariance `C` (typically the converged batch-least-squares
covariance). The state is then parameterized as
`u0 = prior_u0_mean + prior_scale * whitening_L * z`, so that the prior in
whitened coordinates is the standard normal and the posterior is not
ill-conditioned. `prior_u0_std` is ignored in that case.

With `infer_noise = true`, the two noise scales are inferred instead of fixed,
as log-normal parameters with median and log-standard-deviation given by
`sigma_range_prior` and `sigma_doppler_prior`. Fixed scales are appropriate
when the sensor noise is known; inferred scales are appropriate for real data
and for simulation-based calibration of the noise model.
"""
@model function bayesian_orbit_determination_model(
    t_obs, range_obs, doppler_obs, gs::GroundStation,
    prior_u0_mean::Vector{Float64}, prior_u0_std::Vector{Float64},
    opts::OrbitPropagatorOptions,
    whitening_L::Union{Nothing, Matrix{Float64}} = nothing,
    prior_scale::Float64 = 1.0,
    sigma_range::Float64 = 0.010,
    sigma_doppler::Float64 = 1e-4,
    infer_noise::Bool = false,
    sigma_range_prior::Tuple{Float64, Float64} = (0.010, 0.5),
    sigma_doppler_prior::Tuple{Float64, Float64} = (1e-4, 0.5),
)
    z1 ~ Normal(0.0, 1.0)
    z2 ~ Normal(0.0, 1.0)
    z3 ~ Normal(0.0, 1.0)
    z4 ~ Normal(0.0, 1.0)
    z5 ~ Normal(0.0, 1.0)
    z6 ~ Normal(0.0, 1.0)
    z = [z1, z2, z3, z4, z5, z6]

    u0 = if whitening_L === nothing
        prior_u0_mean + prior_u0_std .* z
    else
        prior_u0_mean + prior_scale .* (whitening_L * z)
    end
    u1, u2, u3 = u0[1], u0[2], u0[3]

    # Physical admissibility: reject initial states below the Earth's surface.
    # This is a prior constraint, not part of the force model.
    if norm(SVector(u1, u2, u3)) < R_EARTH + 100.0
        Turing.@addlogprob!(-Inf)
    end

    sigma_r, sigma_d = if infer_noise
        sr ~ LogNormal(log(sigma_range_prior[1]), sigma_range_prior[2])
        sd ~ LogNormal(log(sigma_doppler_prior[1]), sigma_doppler_prior[2])
        (sr, sd)
    else
        (sigma_range, sigma_doppler)
    end

    tspan = (0.0, t_obs[end])
    sol = propagate_orbit(u0, tspan, t_obs, opts)

    for (i, t) in enumerate(t_obs)
        r_sat = SVector{3}(sol.u[i][1], sol.u[i][2], sol.u[i][3])
        v_sat = SVector{3}(sol.u[i][4], sol.u[i][5], sol.u[i][6])
        meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t)

        range_obs[i] ~ Normal(meas.range, sigma_r)
        doppler_obs[i] ~ Normal(meas.range_rate, sigma_d)
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
- `whitening_covariance`: optional 6x6 covariance used to whiten the state
  parameterization (typically the batch-least-squares covariance). Needed when
  the posterior is ill-conditioned, which is the normal case for sparse tracking.
- `prior_scale`: prior width in units of the whitening covariance's standard
  deviations when whitening is used
- `sigma_range`, `sigma_doppler`: standard deviations of the range (km) and
  Doppler (km/s) likelihoods, used when `infer_noise = false`
- `infer_noise`: sample the two noise scales instead of fixing them
- `sigma_range_prior`, `sigma_doppler_prior`: `(median, log_sd)` of the
  log-normal priors on the noise scales; only used when `infer_noise = true`

# Returns
An MCMC chain with samples for `z1, ..., z6`, plus `sr` and `sd` when
`infer_noise = true`. Without whitening, `u0 = prior_mean + prior_std .* z`.
With whitening, `u0 = prior_mean + prior_scale .* cholesky(whitening_covariance).L * z`.
"""
function fit_bayesian_od(
    t_obs::Vector{Float64}, range_obs::Vector{Float64}, doppler_obs::Vector{Float64},
    gs::GroundStation, prior_mean::Vector{Float64}, prior_std::Vector{Float64},
    opts::OrbitPropagatorOptions;
    n_samples::Int = 500, n_adapt::Int = 250, target_accept::Float64 = 0.90,
    whitening_covariance::Union{Nothing, AbstractMatrix{Float64}} = nothing,
    prior_scale::Float64 = 1.0,
    sigma_range::Float64 = 0.010,
    sigma_doppler::Float64 = 1e-4,
    infer_noise::Bool = false,
    sigma_range_prior::Tuple{Float64, Float64} = (0.010, 0.5),
    sigma_doppler_prior::Tuple{Float64, Float64} = (1e-4, 0.5),
)
    length(t_obs) == length(range_obs) == length(doppler_obs) ||
        throw(DimensionMismatch("t_obs, range_obs and doppler_obs must have equal length"))

    whitening_L = if whitening_covariance === nothing
        nothing
    else
        Matrix(cholesky(Symmetric(whitening_covariance)).L)
    end

    model = bayesian_orbit_determination_model(
        t_obs, range_obs, doppler_obs, gs, prior_mean, prior_std, opts,
        whitening_L, prior_scale, sigma_range, sigma_doppler,
        infer_noise, sigma_range_prior, sigma_doppler_prior,
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