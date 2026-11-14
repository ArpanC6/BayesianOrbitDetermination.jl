# Simulation-based calibration of the whitened Bayesian orbit-determination
# posterior, including the two observation-noise scales.
#
# Every number in benchmarks/SBC.md and every mark in
# benchmarks/sbc_rank_histogram.svg is written by this script.
#
# Full run (100 simulations, about 2.5 h on one core):
#     julia --project=. benchmarks/run_sbc_study.jl
#
# Sharded run (8 processes, about 25 min on 8 cores):
#     SBC_N=104 SBC_NSHARDS=8 SBC_SHARD=k julia --project=. benchmarks/run_sbc_study.jl
#     SBC_MERGE=1 julia --project=. benchmarks/run_sbc_study.jl
#
# The prior is fixed before any simulation: it is centered on the batch
# least-squares solution of a pilot data set generated from a reference orbit,
# with covariance 10x the pilot BLS covariance. It does not depend on the data
# of any simulation, which is what makes the rank test valid.

using BayesianOrbitDetermination
using Turing
using Statistics
using Printf
using Random
using LinearAlgebra
using StaticArrays
using Distributions

const N_TOTAL  = parse(Int, get(ENV, "SBC_N", "100"))
const N_SHARDS = parse(Int, get(ENV, "SBC_NSHARDS", "1"))
const SHARD    = parse(Int, get(ENV, "SBC_SHARD", "0"))
const MERGE    = get(ENV, "SBC_MERGE", "0") == "1"
const N_SAMPLES = parse(Int, get(ENV, "SBC_SAMPLES", "500"))
const N_ADAPT   = parse(Int, get(ENV, "SBC_ADAPT", "500"))
const PRIOR_SCALE = parse(Float64, get(ENV, "SBC_PRIOR_SCALE", "10.0"))
const SEED = 2025
const SR_PRIOR = (0.010, 0.5)
const SD_PRIOR = (1e-4, 0.5)

const OUT_DIR = @__DIR__
const SHARD_DIR = joinpath(OUT_DIR, "sbc_shards")

function pilot_prior()
    gs = GroundStation("BenchmarkStation", 35.0, -115.0, 0.5)
    opts = OrbitPropagatorOptions(include_j2 = true, include_drag = true,
                                  Cd = 2.2, area_to_mass = 0.01)
    u_pilot = let c = kepler_to_cartesian(
            KeplerianState(6900.0, 0.002, deg2rad(45.0), 0.3, 1.1, 2.0))
        [c.r...; c.v...]
    end

    t_pilot, n_passes = select_pass_times(u_pilot, gs, opts)
    n_passes >= 2 || error("pilot orbit produced fewer than two passes")

    Random.seed!(123)
    data = simulate_tracking_data(u_pilot, t_pilot, gs, 0.010, 1e-4; opts = opts)
    guess = u_pilot + [0.3, -0.2, 0.1, 1e-4, -2e-4, 1e-4]
    bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                                  guess, 0.010, 1e-4, opts)
    bls.converged || error("pilot BLS did not converge (rms = $(bls.rms_residual))")

    P = Symmetric(Matrix(bls.covariance))
    return (gs = gs, opts = opts, prior_mean = bls.u_estimated,
            prior_covariance = Matrix(P), n_pilot = length(data.t_obs),
            pilot_rms = bls.rms_residual, pilot_cond = cond(Matrix(P)),
            pilot_sigma = sqrt.(diag(P)))
end

function write_shard(res::SBCResult, path::String)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sim,param,rank,truth,median,rhat,ess,n_obs,divergences,seconds")
        for i in 1:size(res.ranks, 1)
            for (j, name) in enumerate(res.param_names)
                println(io, join((i, String(name), res.ranks[i, j],
                                  res.truths[i, j], res.medians[i, j],
                                  res.rhat[i, j], res.ess[i, j],
                                  res.n_obs[i], res.divergences[i], res.seconds[i]), ','))
            end
        end
    end
end

function read_shards()
    files = sort(filter(f -> startswith(basename(f), "rank_") && endswith(f, ".csv"),
                        readdir(SHARD_DIR; join = true)))
    isempty(files) && error("no shard files in $SHARD_DIR")
    # Key by (shard, sim): the CSV row index is local to a shard, so two shards
    # both writing rows 1..n would otherwise overwrite each other.
    rows = Dict{Tuple{Int, Int, String}, Vector{Float64}}()
    for f in files
        shard = parse(Int, match(r"rank_(\d+)\.csv$", basename(f)).captures[1])
        for line in Iterators.drop(eachline(f), 1)
            isempty(strip(line)) && continue
            c = split(line, ',')
            key = (shard, parse(Int, c[1]), c[2])
            rows[key] = parse.(Float64, c[3:end])
        end
    end
    sims = sort(unique((k[1], k[2]) for k in keys(rows)))
    params = sort(unique(k[3] for k in keys(rows)))
    n, d = length(sims), length(params)
    ranks = zeros(Int, n, d)
    truths = zeros(n, d); medians = zeros(n, d)
    rhat = zeros(n, d); ess = zeros(n, d)
    n_obs = zeros(Int, n); div = zeros(Int, n); secs = zeros(n)
    for (i, s) in enumerate(sims), (j, p) in enumerate(params)
        v = rows[(s[1], s[2], p)]
        ranks[i, j] = round(Int, v[1]); truths[i, j] = v[2]; medians[i, j] = v[3]
        rhat[i, j] = v[4]; ess[i, j] = v[5]
        n_obs[i] = round(Int, v[6]); div[i] = round(Int, v[7]); secs[i] = v[8]
    end
    return (sims = sims, params = params, ranks = ranks, truths = truths,
            medians = medians, rhat = rhat, ess = ess, n_obs = n_obs,
            divergences = div, seconds = secs, files = files)
end

function svg_histogram(path, ranks, names, pvals; n_bins = 10, n_draws = 500)
    n, d = size(ranks)
    ncols = 4
    nrows = cld(d, ncols)
    pw, ph = 280, 170
    ml, mt, mb = 46, 34, 40
    W = ncols * pw + 20
    H = nrows * ph + 70
    expected = n / n_bins
    band_lo, band_hi = quantile(Binomial(n, 1 / n_bins), [0.025, 0.975])

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" viewBox="0 0 $W $H" font-family="sans-serif">""")
    println(io, """<rect width="$W" height="$H" fill="white"/>""")
    println(io, """<text x="20" y="26" font-size="16" font-weight="bold">SBC rank histograms (N = $n simulations, $n_draws draws)</text>""")
    println(io, """<text x="20" y="44" font-size="11" fill="#555">Uniform is calibrated. Bars show observed rank counts; grey band is the 95% binomial interval.</text>""")

    for (j, name) in enumerate(names)
        col = (j - 1) % ncols
        row = (j - 1) ÷ ncols
        x0 = 20 + col * pw
        y0 = 64 + row * ph
        counts = zeros(Int, n_bins)
        for r in ranks[:, j]
            counts[min(n_bins, floor(Int, r / (n_draws + 1) * n_bins) + 1)] += 1
        end
        ymax = max(maximum(counts), band_hi) * 1.12
        plot_w, plot_h = pw - ml - 16, ph - mt - mb
        sx = plot_w / n_bins
        sy = plot_h / ymax

        println(io, """<g transform="translate($x0,$y0)">""")
        println(io, """<text x="$ml" y="16" font-size="12" font-weight="bold">$(name)</text>""")
        println(io, """<text x="$(pw - 16)" y="16" font-size="10" text-anchor="end" fill="#555">chi2 p = $(@sprintf("%.3f", pvals[j]))</text>""")
        # axes
        println(io, """<line x1="$ml" y1="$(mt + plot_h)" x2="$(ml + plot_w)" y2="$(mt + plot_h)" stroke="#333"/>""")
        println(io, """<line x1="$ml" y1="$mt" x2="$ml" y2="$(mt + plot_h)" stroke="#333"/>""")
        # 95% band
        by_lo = mt + plot_h - band_lo * sy
        by_hi = mt + plot_h - band_hi * sy
        println(io, """<rect x="$ml" y="$by_hi" width="$plot_w" height="$(by_lo - by_hi)" fill="#e6e6e6"/>""")
        # expected line
        ey = mt + plot_h - expected * sy
        println(io, """<line x1="$ml" y1="$ey" x2="$(ml + plot_w)" y2="$ey" stroke="#111" stroke-dasharray="4,3"/>""")
        # bars
        for b in 1:n_bins
            bx = ml + (b - 1) * sx + 1
            bh = counts[b] * sy
            println(io, """<rect x="$(@sprintf("%.1f", bx))" y="$(@sprintf("%.1f", mt + plot_h - bh))" width="$(@sprintf("%.1f", sx - 2))" height="$(@sprintf("%.1f", bh))" fill="#3b6ea5"/>""")
        end
        # y ticks 0 and expected
        println(io, """<text x="$(ml - 6)" y="$(mt + plot_h + 4)" font-size="9" text-anchor="end" fill="#555">0</text>""")
        println(io, """<text x="$(ml - 6)" y="$(ey + 3)" font-size="9" text-anchor="end" fill="#555">$(@sprintf("%.0f", expected))</text>""")
        # x ticks
        println(io, """<text x="$ml" y="$(mt + plot_h + 15)" font-size="9" fill="#555">0</text>""")
        println(io, """<text x="$(ml + plot_w)" y="$(mt + plot_h + 15)" font-size="9" text-anchor="end" fill="#555">$n_draws</text>""")
        println(io, """<text x="$(ml + plot_w / 2)" y="$(mt + plot_h + 30)" font-size="10" text-anchor="middle" fill="#555">rank of true value</text>""")
        println(io, "</g>")
    end
    println(io, "</svg>")
    open(path, "w") do f
        write(f, String(take!(io)))
    end
end

function merge_and_report(pilot)
    d = read_shards()
    n, np = size(d.ranks)

    # Canonical parameter order
    order = [findfirst(==(String(p)), d.params) for p in
             vcat(["z1", "z2", "z3", "z4", "z5", "z6"], ["log_sigma_range", "log_sigma_doppler"])]
    if any(isnothing, order)
        order = collect(1:np)
    end
    ranks = d.ranks[:, order]
    names = d.params[order]
    ut = sbc_uniformity_test(ranks; n_bins = 10, n_draws = N_SAMPLES)
    truths = d.truths[:, order]
    medians = d.medians[:, order]
    rhat = d.rhat[:, order]
    ess = d.ess[:, order]

    open(joinpath(OUT_DIR, "sbc_ranks.csv"), "w") do io
        println(io, join(vcat(["sim"], names), ','))
        for i in 1:n
            println(io, join(vcat(["$(d.sims[i][1])-$(d.sims[i][2])"], ranks[i, :]), ','))
        end
    end
    svg_histogram(joinpath(OUT_DIR, "sbc_rank_histogram.svg"), ranks, names, ut.p_value;
                  n_bins = 10, n_draws = N_SAMPLES)

    n_div = sum(d.divergences)
    n_div_sims = count(>(0), d.divergences)
    open(joinpath(OUT_DIR, "SBC.md"), "w") do io
        println(io, "# Simulation-Based Calibration")
        println(io)
        println(io, "Generated by `benchmarks/run_sbc_study.jl`. Do not edit by hand.")
        println(io, "Raw ranks: `sbc_ranks.csv`. Figure: `sbc_rank_histogram.svg`.")
        println(io)
        println(io, "## Protocol")
        println(io)
        println(io, "- $(n) simulations, seeded independently as `2025 + i`")
        println(io, "- Reference orbit: a = 6900 km, e = 0.002, i = 45 deg; station 35N, 115W, 0.5 km")
        println(io, "- Force model: two-body + J2 + exponential drag")
        println(io, "- Truth drawn from the fixed prior `N(u_ref, $(PRIOR_SCALE)^2 L L')` with `L L' = P_BLS` of a pilot fit ($(pilot.n_pilot) observations, rms $(@sprintf("%.3f", pilot.pilot_rms)))")
        println(io, "- Noise scales drawn from LogNormal(median $(SR_PRIOR[1]) km, log-sd $(SR_PRIOR[2])) and LogNormal(median $(SD_PRIOR[1]) km/s, log-sd $(SD_PRIOR[2]))")
        println(io, "- Sparse tracking: up to 3 passes over 24 h above 10 deg elevation, 8 range and Doppler samples per pass, sampling times taken from the true orbit")
        println(io, "- Posterior: NUTS, $(N_SAMPLES) draws after $(N_ADAPT) adaptation, target acceptance 0.90, state whitened by `P_BLS`, prior scale $(PRIOR_SCALE)")
        println(io, "- Priors in the fit are exactly the priors used to simulate; the pilot prior is independent of every simulated data set")
        println(io, "- Software: Julia $(VERSION), Turing $(pkgversion(Turing)), BayesianOrbitDetermination $(pkgversion(BayesianOrbitDetermination))")
        println(io)
        println(io, "Pilot prior: condition number $(@sprintf("%.2e", pilot.pilot_cond)); per-axis 1-sigma $(round.(pilot.pilot_sigma, sigdigits = 3)).")
        println(io)
        println(io, "## Rank uniformity")
        println(io)
        println(io, "A calibrated posterior gives uniform ranks on `0:$(N_SAMPLES)` (the number of posterior draws). The chi-square test uses 10 equal-probability rank bins and 9 degrees of freedom.")
        println(io)
        println(io, "| Parameter | chi2 (9 dof) | p-value | median rank / (N+1) | median ESS | max r-hat |")
        println(io, "|---|---|---|---|---|---|")
        for j in 1:length(names)
            @printf(io, "| %s | %.2f | %.3f | %.3f | %.0f | %.3f |\n",
                    names[j], ut.chi2[j], ut.p_value[j],
                    median(ranks[:, j]) / (n + 1), median(ess[:, j]), maximum(rhat[:, j]))
        end
        println(io)
        n_sig = count(p -> p < 0.05, ut.p_value)
        println(io, "$(n_sig) of $(length(names)) parameters have p < 0.05 (expected about $(round(0.05 * length(names), digits = 1)) under calibration); the Bonferroni threshold is $(round(0.05 / length(names), sigdigits = 2)).")
        println(io)
        println(io, "## Sampler diagnostics")
        println(io)
        println(io, "| Item | Value |")
        println(io, "|---|---|")
        println(io, "| Simulations | $(n) |")
        println(io, "| Simulations with divergent transitions | $(n_div_sims) / $(n) |")
        println(io, "| Divergent transitions, total | $(n_div) |")
        @printf(io, "| Worst r-hat | %.3f |\n", maximum(rhat))
        @printf(io, "| Smallest ESS | %.0f |\n", minimum(ess))
        @printf(io, "| Median ESS | %.0f |\n", median(ess))
        @printf(io, "| Observations per simulation | %s |\n", "$(minimum(d.n_obs))-$(maximum(d.n_obs))")
        @printf(io, "| Mean sampling time per simulation | %.0f s |\n", mean(d.seconds))
        println(io, "| Shard files | $(length(d.files)) |")
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "Uniform ranks mean the posterior is calibrated in the prior-to-posterior\n",
                    "sense: neither overconfident nor underconfident, for the state and for the\n",
                    "two noise scales. A U-shaped histogram means the posterior is too narrow;\n",
                    "a dome means it is too wide. The test has limited power at N = $(n): it\n",
                    "detects gross miscalibration, not a few-percent error.")
    end

    @printf("Merged %d shards: %d simulations, %d parameters\n", length(d.files), n, length(names))
    println("Wrote SBC.md, sbc_ranks.csv, sbc_rank_histogram.svg")
end

function main()
    pilot = pilot_prior()
    @printf("Pilot: %d obs, rms %.3f, cond %.2e\n",
            pilot.n_pilot, pilot.pilot_rms, pilot.pilot_cond)

    if MERGE
        merge_and_report(pilot)
        return
    end

    shard = SHARD == 0 ? 1 : SHARD
    nshards = SHARD == 0 ? N_SHARDS : N_SHARDS
    indices = collect(shard:nshards:N_TOTAL)
    println("Shard $(shard)/$(nshards): $(length(indices)) simulations of $(N_TOTAL)")

    t0 = time()
    res = run_sbc(pilot.gs, pilot.opts, pilot.prior_mean, pilot.prior_covariance;
                  n_simulations = N_TOTAL, sim_indices = indices,
                  prior_scale = PRIOR_SCALE,
                  infer_noise = true,
                  sigma_range_prior = SR_PRIOR,
                  sigma_doppler_prior = SD_PRIOR,
                  n_samples = N_SAMPLES, n_adapt = N_ADAPT,
                  target_accept = 0.90, seed = SEED, verbose = true)
    elapsed = time() - t0

    path = joinpath(SHARD_DIR, @sprintf("rank_%02d.csv", shard))
    write_shard(res, path)
    @printf("Shard %d done in %.0f s (%d kept, %d skipped). Wrote %s\n",
            shard, elapsed, size(res.ranks, 1), res.skipped, path)

    if N_SHARDS == 1 && SHARD == 0
        merge_and_report(pilot)
    end
end

main()
