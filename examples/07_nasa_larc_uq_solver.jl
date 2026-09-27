# Tutorial 07: NASA LaRC Uncertainty Quantification (UQ) Challenge Problem Benchmark

using BayesianOrbitDetermination
using Printf
using Random

println("=== BayesianOrbitDetermination.jl Tutorial 07 ===")
println("Target: NASA LaRC UQ Challenge Problem Epistemic vs Aleatory Bounds\n")

u_true = [6900.0, 0.0, 0.0, 0.0, 7.5, 0.0]
n_draws = 500

# Synthetic posterior samples representing mixed uncertainty
samples = randn(n_draws, 6) .* 0.02 .+ u_true'

uq_res = solve_nasa_larc_uq_subproblem(samples, u_true)

@printf("Epistemic 95%% Position Error Bounds: [%.4f, %.4f] km\n", uq_res.epistemic_bounds[1], uq_res.epistemic_bounds[2])
@printf("Aleatory Standard Deviation: %.4f km\n", uq_res.aleatory_std)
@printf("Mission Reliability (Error < 100m): %.2f%%\n", uq_res.reliability_index * 100.0)

println("Tutorial 07 complete!")
