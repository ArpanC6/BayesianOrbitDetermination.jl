module BayesianOrbitDetermination

using LinearAlgebra
using Statistics
using Random
using Distributions
using StaticArrays
using DifferentialEquations
using Turing
using ForwardDiff
using Printf
import SatelliteToolboxSgp4
import SatelliteToolboxTle

# Core Earth Constants (WGS-84 / GRS-80 baseline)
const MU_EARTH = 398600.4418 # km^3 / s^2
const R_EARTH  = 6378.137    # Earth equatorial radius (km)
const J2_EARTH = 1.0826269e-3 # Earth J2 harmonic coefficient
const OMEGA_EARTH = 7.292115e-5 # Earth rotation rate (rad/s)

# Submodules & Files
include("dynamics/two_body.jl")
include("dynamics/perturbations.jl")
include("dynamics/earth_orientation.jl")
include("dynamics/spherical_harmonics.jl")
include("dynamics/sgp4_tle.jl")
include("dynamics/propagator.jl")

include("observations/measurement_models.jl")
include("observations/ingestion.jl")
include("observations/ccsds_oem.jl")
include("observations/isro_stations.jl")

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
export zonal_harmonic_perturbation, j2_perturbation, drag_perturbation, srp_perturbation
export sun_direction_eci, sun_distance_au, gmst_rad, itrf_to_eci_state
export orbit_dynamics, propagate_orbit, OrbitPropagatorOptions
export GroundStation, ground_station_ecef, ecef_to_eci, eci_to_station_azel_range_doppler, simulate_tracking_data, select_pass_times
export GroundStationObservation, parse_ccsds_oem
export CCSDSOEMEphemeris, export_ccsds_oem, parse_ccsds_oem_file
export ISRO_GROUND_STATIONS, get_isro_station
export CPFPrediction, parse_cpf, parse_slr_normal_points, date_to_mjd
export bayesian_orbit_determination_model, fit_bayesian_od
export fit_batch_least_squares
export run_ekf_orbit_determination, run_ukf_orbit_determination, FilterResult, propagate_state_and_stm
export compute_mahalanobis_distance, EvaluateCovarianceRealism, SBCResult, run_sbc, sbc_uniformity_test
export ConjunctionEvent, compute_collision_probability_foster, compute_collision_probability_monte_carlo, export_cdm_summary, export_ccsds_cdm
export NASAUQChallengeResult, solve_nasa_larc_uq_subproblem

end # module BayesianOrbitDetermination
