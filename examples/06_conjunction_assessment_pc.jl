# Tutorial 06: Conjunction Risk Assessment (Pc) and Summary Export
#
# The exported file is a human-readable key-value summary, not a CCSDS CDM.

using BayesianOrbitDetermination
using LinearAlgebra
using Printf
using StaticArrays # Added missing import

println("--- BayesianOrbitDetermination.jl Tutorial 06 ---")
println("Target: Space Collision Probability (Pc) and summary export\n")

# Satellite 1 (ISS)
r1 = SVector(6800.0, 0.0, 0.0)
v1 = SVector(0.0, 7.6, 0.0)
cov1 = Matrix(Diagonal([0.01, 0.01, 0.01, 1e-6, 1e-6, 1e-6]))

# Satellite 2 (Debris object passing within 150 meters)
r2 = SVector(6800.12, 0.05, 0.02)
v2 = SVector(0.0, -7.6, 1.2) # High relative velocity encounter
cov2 = Matrix(Diagonal([0.05, 0.05, 0.05, 5e-6, 5e-6, 5e-6]))

conj = ConjunctionEvent("ISS", "COSMOS-2251-DEBRIS", 1000.0, r1, v1, cov1, r2, v2, cov2, 20.0) # 20m combined hard-body radius

pc = compute_collision_probability_foster(conj)
@printf("Encounter Miss Distance: %.3f km (%.1f meters)\n", norm(r2 - r1), norm(r2 - r1)*1000.0)
@printf("2D Foster-Elrod Collision Probability Pc: %.6e\n", pc)

summary_file = joinpath(@__DIR__, "sample_conjunction_summary.txt")
export_cdm_summary(conj, pc, summary_file)
println("Generated conjunction summary (not a CCSDS CDM): $(summary_file)")
println("Tutorial 06 complete")