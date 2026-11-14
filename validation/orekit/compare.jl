# Compare the Julia propagator against the Orekit-generated reference states.
# Run from the repository root:
#   julia --project=. validation/orekit/compare.jl

using BayesianOrbitDetermination
using DelimitedFiles
using LinearAlgebra
using Printf
using Statistics

const U0 = [3205.830222953533, 4997.609661816432, 3488.022461182947,
            -6.035351552351067, 0.6629103392316869, 4.597259758177123]

function compare_case(name, path; drag::Bool)
    data = readdlm(path, ',', Float64; skipstart = 2)
    t = data[:, 1]
    ref = data[:, 2:7]

    opts = OrbitPropagatorOptions(include_j2 = true, include_drag = drag,
                                  reltol = 1e-12, abstol = 1e-12)
    sol = propagate_orbit(U0, (0.0, t[end]), collect(t), opts)
    got = reduce(vcat, permutedims.(collect.(sol.u)))

    dp = [norm(got[i, 1:3] .- ref[i, 1:3]) for i in eachindex(t)]
    dv = [norm(got[i, 4:6] .- ref[i, 4:6]) for i in eachindex(t)]

    @printf("%-7s %2d points up to %6.0f s | max |dr| = %8.4f m | final |dr| = %8.4f m | max |dv| = %10.6f mm/s\n",
            name, length(t), t[end], 1000 * maximum(dp), 1000 * dp[end], 1e6 * maximum(dv))
    return maximum(dp), dp[end], maximum(dv)
end

println("Orekit 13.1.9 cross-validation (EME2000, Dormand-Prince 8(5,3), absTol=1e-6 m)")
println("Julia: Tsit5, reltol = abstol = 1e-12")
compare_case("J2", joinpath(@__DIR__, "reference", "reference_j2.csv"); drag = false)
compare_case("J2+drag", joinpath(@__DIR__, "reference", "reference_j2drag.csv"); drag = true)
