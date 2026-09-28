# Tutorial 06: Conjunction Risk Assessment (Pc) & CCSDS CDM Export

using BayesianOrbitDetermination
using LinearAlgebra
using Printf
using StaticArrays # Added missing import

println("=== BayesianOrbitDetermination.jl Tutorial 06 ===")
println("Target: Space Collision Probability (Pc) & CCSDS CDM Generation\n")

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

cdm_file = joinpath(@__DIR__, "sample_conjunction.cdm")
export_ccsds_cdm(conj, pc, cdm_file)
println("Generated CCSDS Conjunction Data Message: $(cdm_file)")
println("Tutorial 06 complete")