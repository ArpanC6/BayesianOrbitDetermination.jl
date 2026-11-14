# Simulation-Based Calibration (SBC) for orbit determination.
#
# If the prior and likelihood are correct and the posterior is sampled exactly,
# the rank of the true parameter among posterior draws is uniform on
# {0, ..., n} (Talts et al., 2018). The prior used here must be fixed and
# independent of the simulated data; centering it on a fit to the same data
# invalidates the calibration argument.

const SBC_STATE_NAMES = (:z1, :z2, :z3, :z4, :z5, :z6)

"""
    SBCResult

Output of [`run_sbc`](@ref). Fields are indexed by simulation (rows) and
parameter (columns).

- `ranks`       : rank of the true value among posterior draws, `0:n_samples`
- `param_names` : parameter names, `z1..z6` and, with noise inference,
  `log_sigma_range`, `log_sigma_doppler`
- `truths`      : simulated true values in the sampled parameterization
- `medians`     : posterior medians
- `rhat`, `ess` : per-parameter convergence diagnostics
- `divergences` : divergent transitions per simulation
- `n_obs`       : number of observations per simulation
- `seconds`     : wall time per simulation
- `n_draws`     : posterior draws per simulation (ranks lie in `0:n_draws`)
- `skipped`     : simulations dropped before sampling (Earth intersection,
  too few observations)
"""
struct SBCResult
    ranks::Matrix{Int}
    param_names::Vector{Symbol}
    truths::Matrix{Float64}
    medians::Matrix{Float64}
    rhat::Matrix{Float64}
    ess::Matrix{Float64}
    divergences::Vector{Int}
    n_obs::Vector{Int}
    seconds::Vector{Float64}
    n_draws::Int
    skipped::Int
end

"""
    run_sbc(gs, opts, prior_mean, prior_covariance;
            n_simulations = 100, sim_indices = 1:n_simulations,
            prior_scale = 10.0, infer_noise = true,
            sigma_range_prior = (0.010, 0.5),
            sigma_doppler_prior = (1e-4, 0.5),
            n_samples = 500, n_adapt = 500, target_accept = 0.90,
            t_end = 24 * 3600.0, dt = 10.0, per_pass = 8, max_passes = 3,
            min_elevation_deg = 10.0, seed = 2025, verbose = true)

Run simulation-based calibration for the whitened orbit-determination
posterior over the initial state and, optionally, the two observation-noise
scales.

For simulation `i`, `z_true ~ Normal(0, I)` is drawn from the prior
`u_true = prior_mean + prior_scale * L * z_true` with `L L' = prior_covariance`;
the noise scales are drawn from their log-normal priors; synthetic sparse
tracking data are generated from `u_true`; and [`fit_bayesian_od`](@ref) is run
with exactly that prior. Ranks of `z_true`, `log(sigma_range_true)` and
`log(sigma_doppler_true)` among the posterior draws are returned.

`prior_mean` and `prior_covariance` must be fixed before any data are
simulated. A prior estimated from the same observations being fitted does not
produce calibrated ranks.

Simulations are seeded independently as `seed + i`, so any subset
`sim_indices` reproduces the same draws and the study can be sharded across
processes.
"""
function run_sbc(
    gs::GroundStation,
    opts::OrbitPropagatorOptions,
    prior_mean::Vector{Float64},
    prior_covariance::AbstractMatrix{Float64};
    n_simulations::Int = 100,
    sim_indices::AbstractVector{Int} = collect(1:n_simulations),
    prior_scale::Float64 = 10.0,
    infer_noise::Bool = true,
    fixed_sigma_range::Float64 = 0.010,
    fixed_sigma_doppler::Float64 = 1e-4,
    sigma_range_prior::Tuple{Float64, Float64} = (0.010, 0.5),
    sigma_doppler_prior::Tuple{Float64, Float64} = (1e-4, 0.5),
    n_samples::Int = 500,
    n_adapt::Int = 500,
    target_accept::Float64 = 0.90,
    t_end::Float64 = 24 * 3600.0,
    dt::Float64 = 10.0,
    per_pass::Int = 8,
    max_passes::Int = 3,
    min_elevation_deg::Float64 = 10.0,
    min_obs::Int = 6,
    seed::Int = 2025,
    verbose::Bool = true,
)
    length(prior_mean) == 6 || throw(DimensionMismatch("prior_mean must have length 6"))
    size(prior_covariance) == (6, 6) || throw(DimensionMismatch("prior_covariance must be 6x6"))
    all(i -> 1 <= i <= n_simulations, sim_indices) ||
        throw(ArgumentError("sim_indices must lie in 1:n_simulations"))

    L = Matrix(cholesky(Symmetric(Matrix(prior_covariance))).L)
    param_names = infer_noise ?
        vcat(collect(SBC_STATE_NAMES), [:log_sigma_range, :log_sigma_doppler]) :
        collect(SBC_STATE_NAMES)
    diag_names = infer_noise ?
        vcat(collect(SBC_STATE_NAMES), [:sr, :sd]) :
        collect(SBC_STATE_NAMES)
    dim = length(param_names)

    ranks = Int[]
    truths = Float64[]
    medians = Float64[]
    rhats = Float64[]
    esses = Float64[]
    divergences = Int[]
    n_obs = Int[]
    seconds = Float64[]
    skipped = 0

    for i in sim_indices
        Random.seed!(seed + i)

        z_true = randn(6)
        u_true = prior_mean + prior_scale .* (L * z_true)
        if norm(u_true[1:3]) < R_EARTH + 100.0
            skipped += 1
            continue
        end

        log_sr_true = log(sigma_range_prior[1]) + sigma_range_prior[2] * randn()
        log_sd_true = log(sigma_doppler_prior[1]) + sigma_doppler_prior[2] * randn()
        sr_true = exp(log_sr_true)
        sd_true = exp(log_sd_true)

        t_sim, _ = select_pass_times(u_true, gs, opts;
                                     t_end = t_end, dt = dt,
                                     min_elevation_deg = min_elevation_deg,
                                     per_pass = per_pass, max_passes = max_passes)
        if length(t_sim) < min_obs
            skipped += 1
            continue
        end

        t0 = time()
        data = simulate_tracking_data(u_true, t_sim, gs, sr_true, sd_true; opts = opts)
        chain = fit_bayesian_od(
            data.t_obs, data.range, data.doppler, gs,
            prior_mean, sqrt.(diag(Symmetric(prior_covariance))), opts;
            n_samples = n_samples, n_adapt = n_adapt,
            target_accept = target_accept,
            whitening_covariance = Matrix(prior_covariance),
            prior_scale = prior_scale,
            sigma_range = fixed_sigma_range,
            sigma_doppler = fixed_sigma_doppler,
            infer_noise = infer_noise,
            sigma_range_prior = sigma_range_prior,
            sigma_doppler_prior = sigma_doppler_prior,
        )
        elapsed = time() - t0

        z_draws = hcat([vec(chain[k]) for k in SBC_STATE_NAMES]...)
        rank_row = Int[]
        truth_row = Float64[]
        median_row = Float64[]
        for j in 1:6
            push!(rank_row, count(x -> x < z_true[j], z_draws[:, j]))
            push!(truth_row, z_true[j])
            push!(median_row, median(z_draws[:, j]))
        end
        if infer_noise
            for (draws, truth) in ((vec(chain[:sr]), sr_true), (vec(chain[:sd]), sd_true))
                log_draws = log.(draws)
                push!(rank_row, count(x -> x < log(truth), log_draws))
                push!(truth_row, log(truth))
                push!(median_row, log(median(draws)))
            end
        end

        rhat_row = [Turing.rhat(chain)[k] for k in diag_names]
        ess_row = [Turing.ess(chain)[k] for k in diag_names]
        append!(ranks, rank_row)
        append!(truths, truth_row)
        append!(medians, median_row)
        append!(rhats, rhat_row)
        append!(esses, ess_row)
        push!(divergences,
              round(Int, sum(skipmissing(vec(chain[Turing.FlexiChains.Extra(:numerical_error)])))))
        push!(n_obs, length(data.t_obs))
        push!(seconds, elapsed)

        if verbose
            @printf("  sim %3d  n_obs=%2d  rhat=%.3f  ess=%4.0f  div=%d  %.0f s\n",
                    i, length(data.t_obs), maximum(rhat_row),
                    minimum(ess_row), divergences[end], elapsed)
        end
    end

    n_kept = length(seconds)
    return SBCResult(
        permutedims(reshape(ranks, dim, n_kept)),
        param_names,
        permutedims(reshape(truths, dim, n_kept)),
        permutedims(reshape(medians, dim, n_kept)),
        permutedims(reshape(rhats, dim, n_kept)),
        permutedims(reshape(esses, dim, n_kept)),
        divergences, n_obs, seconds, n_samples, skipped,
    )
end

"""
    sbc_uniformity_test(ranks; n_bins = 10, n_draws = maximum(ranks))

Chi-square test of rank uniformity for each parameter. Under a calibrated
posterior the rank of the true value among `n_draws` posterior draws is
uniform on `0:n_draws`, so the bins must be scaled by `n_draws`, not by the
number of simulations. Returns a NamedTuple with per-parameter `chi2`,
`p_value`, plus `dof` and `n_simulations`.
"""
function sbc_uniformity_test(ranks::AbstractMatrix{<:Integer};
                             n_bins::Int = 10,
                             n_draws::Int = max(1, maximum(ranks)))
    n, dim = size(ranks)
    n > 0 || throw(ArgumentError("ranks is empty"))
    n_bins >= 2 || throw(ArgumentError("n_bins must be at least 2"))
    n_draws >= 1 || throw(ArgumentError("n_draws must be positive"))

    chi2 = zeros(dim)
    p_value = zeros(dim)
    expected = n / n_bins
    for j in 1:dim
        counts = zeros(Int, n_bins)
        for r in ranks[:, j]
            counts[min(n_bins, floor(Int, r / (n_draws + 1) * n_bins) + 1)] += 1
        end
        chi2[j] = sum((counts .- expected) .^ 2 ./ expected)
        p_value[j] = ccdf(Chisq(n_bins - 1), chi2[j])
    end
    return (chi2 = chi2, dof = n_bins - 1, p_value = p_value, n_simulations = n)
end
