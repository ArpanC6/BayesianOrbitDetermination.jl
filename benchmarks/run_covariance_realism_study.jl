# Covariance realism study: frequentist BLS vs Bayesian NUTS on identical data.
# Every number in benchmarks/RESULTS.md is written by this script; nothing is typed by hand.

using BayesianOrbitDetermination
using Turing
using Statistics
using Printf
using Random
using LinearAlgebra
using StaticArrays

const N_TRIALS  = parse(Int, get(ENV, "BOD_TRIALS", "50"))
const N_SAMPLES = parse(Int, get(ENV, "BOD_SAMPLES", "500"))
const N_ADAPT   = parse(Int, get(ENV, "BOD_ADAPT", "500"))
const Z90       = 1.6448536269514722
const Z_NAMES   = [:z1, :z2, :z3, :z4, :z5, :z6]

# Sparse tracking means a few short passes, not a continuous arc. Find the
# visibility windows over a day and sample a fixed number of points per pass.
function pass_times(u0, gs, opts; t_end = 24 * 3600.0, dt = 10.0,
                    min_elevation_deg = 10.0, per_pass = 8, max_passes = 3)
    grid = collect(0.0:dt:t_end)
    sol = propagate_orbit(u0, (0.0, t_end), grid, opts)
    min_el = deg2rad(min_elevation_deg)

    visible = falses(length(grid))
    for i in eachindex(grid)
        r_sat = SVector{3}(sol.u[i][1], sol.u[i][2], sol.u[i][3])
        v_sat = SVector{3}(sol.u[i][4], sol.u[i][5], sol.u[i][6])
        visible[i] = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, grid[i]).elevation > min_el
    end

    times = Float64[]
    n_passes = 0
    i = 1
    while i <= length(grid) && n_passes < max_passes
        if !visible[i]
            i += 1
            continue
        end
        j = i
        while j < length(grid) && visible[j + 1]
            j += 1
        end
        idx = unique(round.(Int, range(i, j, length = min(per_pass, j - i + 1))))
        append!(times, grid[idx])
        n_passes += 1
        i = j + 1
    end
    return times, n_passes
end

function main()
    Random.seed!(42)

    gs        = GroundStation("BenchmarkStation", 35.0, -115.0, 0.5)
    opts      = OrbitPropagatorOptions(include_j2 = true, include_drag = true,
                                       Cd = 2.2, area_to_mass = 0.01)
    # Bayesian state is whitened by the BLS covariance: u = u_bls + 10*L*z, so the
    # prior is 10x wider than the BLS covariance and the likelihood dominates.
    # Without whitening the posterior condition number is ~1e11 and NUTS cannot mix.
    prior_scale = 10.0

    d2_bls    = Float64[]
    d2_bayes  = Float64[]
    mae_bls   = Float64[]
    mae_bayes = Float64[]
    cov_bls   = zeros(Int, 6)
    cov_bayes = zeros(Int, 6)

    n_used            = 0
    n_too_few_obs     = 0
    n_bls_unconverged = 0
    n_high_rms        = 0

    n_divergent_chains = 0
    n_divergent_total  = 0
    rhat_worst         = 0.0
    ess_worst          = Inf
    trial_rows         = NamedTuple[]

    println("Running $(N_TRIALS) replications...")

    for i in 1:N_TRIALS
        kepler_true = KeplerianState(6900.0 + randn() * 10.0, 0.002,
                                     deg2rad(rand() * 30.0 + 30.0),
                                     rand() * 2pi, rand() * 2pi, rand() * 2pi)
        cart_true = kepler_to_cartesian(kepler_true)
        u0_true   = [cart_true.r...; cart_true.v...]

        t_pass, n_passes = pass_times(u0_true, gs, opts)
        if n_passes < 2
            n_too_few_obs += 1
            continue
        end

        data = simulate_tracking_data(u0_true, t_pass, gs, 0.010, 1e-4; opts = opts)
        if length(data.t_obs) < 10
            n_too_few_obs += 1
            continue
        end

        guess = u0_true + randn(6) .* [0.5, 0.5, 0.5, 5e-4, 5e-4, 5e-4]

        bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                                      guess, 0.010, 1e-4, opts)
        # Only stationary fits enter the calibration statistics; a stopped
        # optimizer reports a covariance around a point it has not found.
        if !bls.converged
            n_bls_unconverged += 1
            continue
        end
        # A fit is usable when the normalized post-fit residuals are order unity,
        # which is what "the model and the stated noise are consistent" means.
        if !isfinite(bls.rms_residual) || bls.rms_residual > 2.0
            n_high_rms += 1
            continue
        end

        prior_mean = bls.u_estimated
        cov_bls_sym = Symmetric(Matrix(bls.covariance))
        whitening_L = Matrix(cholesky(cov_bls_sym).L)
        chain = fit_bayesian_od(data.t_obs, data.range, data.doppler, gs,
                                prior_mean, sqrt.(diag(cov_bls_sym)), opts;
                                n_samples = N_SAMPLES, n_adapt = N_ADAPT,
                                target_accept = 0.90,
                                whitening_covariance = Matrix(cov_bls_sym),
                                prior_scale = prior_scale)
        z = hcat([vec(chain[k]) for k in Z_NAMES]...)
        u = (prior_mean .+ prior_scale .* (whitening_L * z'))'

        # Divergence flags live in FlexiChains metadata (Turing 0.49+).
        n_div = round(Int, sum(skipmissing(vec(chain[Turing.FlexiChains.Extra(:numerical_error)]))))
        rhat_trial = maximum(Turing.rhat(chain)[k] for k in Z_NAMES)
        ess_trial  = minimum(Turing.ess(chain)[k] for k in Z_NAMES)
        n_divergent_total += n_div
        n_divergent_chains += (n_div > 0)
        rhat_worst = max(rhat_worst, rhat_trial)
        ess_worst  = min(ess_worst, ess_trial)

        ev = EvaluateCovarianceRealism(u0_true, u, cov(u))

        d2_bls_trial   = compute_mahalanobis_distance(u0_true, bls.u_estimated, bls.covariance)
        mae_bls_trial  = mean(abs.(u0_true[1:3] .- bls.u_estimated[1:3]))
        push!(d2_bls,    d2_bls_trial)
        push!(d2_bayes,  ev.mahalanobis_d2)
        push!(mae_bls,   mae_bls_trial)
        push!(mae_bayes, ev.mae_position)
        push!(trial_rows, (trial = i, n_obs = length(data.t_obs), n_passes = n_passes,
                           bls_rms = bls.rms_residual, bls_converged = bls.converged,
                           d2_bls = d2_bls_trial, d2_bayes = ev.mahalanobis_d2,
                           mae_bls_km = mae_bls_trial, mae_bayes_km = ev.mae_position,
                           divergences = n_div, max_rhat = rhat_trial, min_ess = ess_trial))

        for j in 1:6
            s = sqrt(bls.covariance[j, j])
            lo = bls.u_estimated[j] - Z90 * s
            hi = bls.u_estimated[j] + Z90 * s
            cov_bls[j] += (u0_true[j] >= lo && u0_true[j] <= hi)
        end
        cov_bayes .+= ev.coverage_90

        n_used += 1
        if i % 10 == 0
            @printf("  %d / %d trials, %d used\n", i, N_TRIALS, n_used)
        end
    end

    if n_used == 0
        println("No trial produced a usable fit. Nothing written.")
        return
    end

    mean_se(x) = (mean(x), std(x) / sqrt(length(x)))

    d2_b_m,  d2_b_se  = mean_se(d2_bayes)
    d2_m_m,  d2_m_se  = mean_se(d2_bls)
    mae_b_m, mae_b_se = mean_se(mae_bayes)
    mae_m_m, mae_m_se = mean_se(mae_bls)

    bls_cov_pct   = (cov_bls   ./ n_used) .* 100.0
    bayes_cov_pct = (cov_bayes ./ n_used) .* 100.0

    println()
    @printf("Used trials            : %d / %d\n", n_used, N_TRIALS)
    @printf("Skipped (too few obs)  : %d\n", n_too_few_obs)
    @printf("Skipped (BLS not converged): %d\n", n_bls_unconverged)
    @printf("Skipped (BLS rms > 2)   : %d\n", n_high_rms)
    @printf("BLS   D^2              : %.3f +/- %.3f  (ideal 6.0)\n", d2_m_m, d2_m_se)
    @printf("Bayes D^2              : %.3f +/- %.3f  (ideal 6.0)\n", d2_b_m, d2_b_se)
    @printf("BLS   position MAE     : %.4f +/- %.4f km\n", mae_m_m, mae_m_se)
    @printf("Bayes position MAE     : %.4f +/- %.4f km\n", mae_b_m, mae_b_se)
    println("BLS   per-component 90% coverage (%): ", round.(bls_cov_pct, digits = 1))
    println("Bayes per-component 90% coverage (%): ", round.(bayes_cov_pct, digits = 1))
    @printf("Chains with divergences: %d / %d (total %d)\n",
            n_divergent_chains, n_used, n_divergent_total)
    @printf("Worst r-hat            : %.3f\n", rhat_worst)
    @printf("Smallest ESS           : %.0f\n", ess_worst)

    out = get(ENV, "BOD_OUT", joinpath(@__DIR__, "RESULTS.md"))
    csv = joinpath(dirname(abspath(out)), "results_per_trial.csv")
    open(csv, "w") do io
        cols = propertynames(trial_rows[1])
        println(io, join(cols, ","))
        for row in trial_rows
            println(io, join((getproperty(row, c) for c in cols), ","))
        end
    end

    open(out, "w") do io
        println(io, "# Covariance Realism Study")
        println(io)
        println(io, "Generated by `benchmarks/run_covariance_realism_study.jl`. Do not edit by hand.")
        println(io, "Per-trial values are in `results_per_trial.csv`.")
        println(io)
        println(io, "## Protocol")
        println(io)
        println(io, "- LEO orbit, two-body + J2 + exponential drag")
        println(io, "- Sparse tracking: up to 3 passes over 24 h above 10 deg elevation, 8 range (10 m) and Doppler (1 mm/s) samples per pass")
        println(io, "- Truth: a = 6900 +/- 10 km, e = 0.002, i uniform in 30-60 deg, node/argp/anomaly uniform")
        println(io, "- Initial guess perturbed by 0.5 km and 0.5 m/s per axis")
        println(io, "- Bayesian state is whitened by the BLS covariance: u = u_BLS + $(prior_scale) * L * z, so the prior is $(prior_scale)x wider than the BLS covariance")
        println(io, "- $(N_TRIALS) replications, common random seed 42")
        println(io, "- BLS: Levenberg-Marquardt, Gauss-Newton covariance, trials dropped when the fit does not converge or the normalized post-fit residual RMS exceeds 2")
        println(io, "- Bayes: NUTS, $(N_SAMPLES) draws after $(N_ADAPT) adaptation, target acceptance 0.90")
        println(io, "- Software: Julia $(VERSION), Turing $(pkgversion(Turing)), BayesianOrbitDetermination $(pkgversion(BayesianOrbitDetermination))")
        println(io)
        println(io, "## Trials")
        println(io)
        println(io, "| Item | Count |")
        println(io, "|---|---|")
        println(io, "| Used | $(n_used) |")
        println(io, "| Skipped, fewer than 10 observations | $(n_too_few_obs) |")
        println(io, "| Skipped, BLS did not converge | $(n_bls_unconverged) |")
        println(io, "| Skipped, BLS post-fit residual RMS > 2 | $(n_high_rms) |")
        println(io)
        println(io, "## Results")
        println(io)
        println(io, "| Metric | BLS | Bayesian NUTS | Ideal |")
        println(io, "|---|---|---|---|")
        @printf(io, "| Mean Mahalanobis D^2 | %.3f +/- %.3f | %.3f +/- %.3f | 6.000 |\n",
                d2_m_m, d2_m_se, d2_b_m, d2_b_se)
        @printf(io, "| Position MAE (km) | %.4f +/- %.4f | %.4f +/- %.4f | 0 |\n",
                mae_m_m, mae_m_se, mae_b_m, mae_b_se)
        println(io, "| Per-component 90% coverage | $(round.(bls_cov_pct, digits = 1)) | $(round.(bayes_cov_pct, digits = 1)) | 90 |")
        println(io)
        println(io, "## MCMC diagnostics")
        println(io)
        println(io, "| Item | Value |")
        println(io, "|---|---|")
        println(io, "| Chains with divergent transitions | $(n_divergent_chains) / $(n_used) |")
        println(io, "| Divergent transitions, total | $(n_divergent_total) |")
        @printf(io, "| Worst r-hat over z1..z6 | %.3f |\n", rhat_worst)
        @printf(io, "| Smallest bulk ESS over z1..z6 | %.0f |\n", ess_worst)
        println(io)
        println(io, "Uncertainties are one standard error over replications.")
    end

    println("\nWrote ", out)
end

main()
