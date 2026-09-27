module BayesianOrbitDetermination

using LinearAlgebra
using Statistics
using Random
using Distributions
using StaticArrays
using DifferentialEquations
using Turing
using Optim
using ForwardDiff
using Printf

# Core Earth Constants (WGS-84 / GRS-80 baseline)
const MU_EARTH = 398600.4418 # km^3 / s^2
const R_EARTH  = 6378.137    # Earth equatorial radius (km)
const J2_EARTH = 1.0826269e-3 # Earth J2 harmonic coefficient
const OMEGA_EARTH = 7.292115e-5 # Earth rotation rate (rad/s)

# Submodules & Files
include("dynamics/two_body.jl")
include("dynamics/perturbations.jl")
include("dynamics/spherical_harmonics.jl")
include("dynamics/sgp4_tle.jl")
include("dynamics/propagator.jl")

include("observations/measurement_models.jl")
include("observations/ingestion.jl")

include("inference/bayesian_od.jl")
include("inference/optimizers.jl")
include("inference/sequential_filters.jl")

include("realism/covariance_realism.jl")
include("realism/sbc.jl")

include("conjunction/collision_probability.jl")
include("conjunction/cdm_exporter.jl")

include("agency/nasa_larc_uq.jl")

# Exports
export MU_EARTH, R_EARTH, J2_EARTH, OMEGA_EARTH, J3_EARTH, J4_EARTH, J5_EARTH
export KeplerianState, CartesianState, kepler_to_cartesian, cartesian_to_kepler
export TLE, parse_tle, tle_to_keplerian, propagate_sgp4
export zonal_harmonic_perturbation
export orbit_dynamics, propagate_orbit, OrbitPropagatorOptions
export GroundStation, ecef_to_eci, eci_to_station_azel_range_doppler, simulate_tracking_data
export GroundStationObservation, parse_ccsds_oem
export bayesian_orbit_determination_model, fit_bayesian_od
export fit_batch_least_squares
export run_ekf_orbit_determination, run_ukf_orbit_determination, FilterResult
export compute_mahalanobis_distance, EvaluateCovarianceRealism, compute_sbc_rank_statistics
export ConjunctionEvent, compute_collision_probability_foster, export_ccsds_cdm
export NASAUQChallengeResult, solve_nasa_larc_uq_subproblem

end # module BayesianOrbitDetermination
