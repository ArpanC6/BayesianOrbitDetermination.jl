# Example 09: Indian Space Situational Awareness (ISRO NETRA) Simulation
#
# Simulates:
# 1. Indian Earth Observation Satellite (e.g. CARTOSAT-2 / EOS series) in Sun-Synchronous LEO.
# 2. Tracking passes over ISRO ground stations (SDSC SHAR Sriharikota & ISTRAC Bengaluru).
# 3. Covariance-whitened Bayesian Orbit Determination under sparse tracking.
# 4. Standard CCSDS 502.0-B-2 OEM ephemeris generation.
# 5. Conjunction risk analysis against tracked orbital debris: Foster-Elrod Pc,
#    independent Monte Carlo validation, and CCSDS 508.0-B-1 CDM generation.

using BayesianOrbitDetermination
using LinearAlgebra
using StaticArrays
using Printf

println("================================================================================")
println(" ISRO SSA / NETRA TRACKING & CONJUNCTION RISK ASSESSMENT DEMONSTRATION")
println("================================================================================")

# 1. Define Indian Satellite & ISRO Ground Stations
# CARTOSAT-class SSO orbit: alt ~630 km (a = 7008 km), e = 0.001, inc = 97.8 deg
k_cartosat = KeplerianState(7008.137, 0.0014, deg2rad(97.8), deg2rad(45.0), deg2rad(120.0), deg2rad(0.0))
c_cartosat = kepler_to_cartesian(k_cartosat)
u0_true = [c_cartosat.r...; c_cartosat.v...]

shar_station = get_isro_station(:SRIHARIKOTA_SHAR)
bengaluru_station = get_isro_station(:BANGALORE_ISTRAC)

println("[1] Asset: CARTOSAT-2 / EOS SSO Satellite")
@printf("    Initial Position ECI: [%.3f, %.3f, %.3f] km (altitude: %.1f km)\n",
        u0_true[1], u0_true[2], u0_true[3], norm(u0_true[1:3]) - R_EARTH)
println("    Ground Stations:")
@printf("    - %s: lat = %.4f°, lon = %.4f°, alt = %.3f km\n",
        shar_station.name, shar_station.latitude_deg, shar_station.longitude_deg, shar_station.altitude_km)
@printf("    - %s: lat = %.4f°, lon = %.4f°, alt = %.3f km\n",
        bengaluru_station.name, bengaluru_station.latitude_deg, bengaluru_station.longitude_deg, bengaluru_station.altitude_km)

# 2. Simulate Radar Tracking Passes
opts = OrbitPropagatorOptions(include_j2 = true, Cd = 2.2, area_to_mass = 0.012)
t_pass_shar, n_passes_shar = select_pass_times(u0_true, shar_station, opts;
                                              t_end = 12 * 3600.0, dt = 15.0, per_pass = 8, max_passes = 2)

@printf("\n[2] Tracking Simulation:\n")
@printf("    Identified %d pass windows over Sriharikota (%d tracking points)\n", n_passes_shar, length(t_pass_shar))

sigma_range_radar = 0.015    # 15 m radar range accuracy
sigma_doppler_radar = 2.0e-4 # 0.2 m/s Doppler accuracy
meas = simulate_tracking_data(u0_true, t_pass_shar, shar_station, sigma_range_radar, sigma_doppler_radar; opts = opts)

# 3. Fit Batch Least Squares
guess = u0_true + [0.4, -0.3, 0.2, 4e-4, -3e-4, 2e-4]
bls = fit_batch_least_squares(meas.t_obs, meas.range, meas.doppler, shar_station,
                              guess, sigma_range_radar, sigma_doppler_radar, opts)

pos_err_bls = norm(bls.u_estimated[1:3] - u0_true[1:3])
@printf("\n[3] Levenberg-Marquardt Batch Least Squares Estimation:\n")
@printf("    Converged: %s | Iterations: %d | Post-fit RMS: %.3f\n", bls.converged, bls.iterations, bls.rms_residual)
@printf("    Position Error from Truth: %.2f meters\n", pos_err_bls * 1000.0)

# 4. Export CCSDS OEM Ephemeris
println("\n[4] Generating CCSDS 502.0-B-2 Orbit Ephemeris Message (OEM)...")
t_grid = collect(0.0:60.0:3600.0)
sol_ephem = propagate_orbit(bls.u_estimated, (0.0, 3600.0), t_grid, opts)
states_mat = zeros(length(t_grid), 6)
epochs_str = String[]
for (i, t) in enumerate(t_grid)
    states_mat[i, :] = sol_ephem.u[i]
    push!(epochs_str, @sprintf("2026-10-29T12:%02d:%02.1f", div(Int(t), 60), mod(t, 60.0)))
end

oem = CCSDSOEMEphemeris("CARTOSAT-2", "2007-001A", "EARTH", "EME2000", "UTC",
                        epochs_str, t_grid, states_mat)
oem_file = joinpath(@__DIR__, "cartosat_ephemeris.oem")
export_ccsds_oem(oem_file, oem)
println("    ✓ Exported standard CCSDS OEM: $oem_file")

# 5. Conjunction Risk Assessment against Debris
println("\n[5] Conjunction Assessment & Collision Probability (Pc):")
r_debris = c_cartosat.r + SVector(0.040, 0.020, -0.030) # ~53 m miss distance
v_debris = c_cartosat.v + SVector(0.5, -14.2, 0.3)      # Hyperbolic relative crossing ~14 km/s
cov_cartosat = Matrix(bls.covariance)
cov_debris = Matrix(Diagonal([0.05^2, 0.05^2, 0.05^2, 1e-4^2, 1e-4^2, 1e-4^2]))

conj_event = ConjunctionEvent("CARTOSAT-2", "COSMOS-DEBRIS", 3600.0,
                              c_cartosat.r, c_cartosat.v, cov_cartosat,
                              r_debris, v_debris, cov_debris, 12.0) # 12 m combined hard-body sphere

pc_foster = compute_collision_probability_foster(conj_event)
pc_mc, pc_mc_se = compute_collision_probability_monte_carlo(conj_event; n_samples = 50_000)

@printf("    Miss Distance: %.2f meters\n", norm(r_debris - c_cartosat.r) * 1000.0)
@printf("    Relative Speed: %.2f km/s\n", norm(v_debris - c_cartosat.v))
@printf("    Foster-Elrod Pc (2D Quad): %.4e\n", pc_foster)
@printf("    Monte Carlo Pc (50k draws): %.4e +/- %.4e\n", pc_mc, pc_mc_se)

# 6. Export Conjunction Data Message (CCSDS 508.0-B-1 CDM)
cdm_file = joinpath(@__DIR__, "cartosat_conjunction.cdm")
export_ccsds_cdm(conj_event, pc_foster, cdm_file;
                 tca_calendar = "2026-10-29T13:00:00.000",
                 message_id = "ISRO-IS4OM-20261029-001",
                 originator = "IS4OM / BayesianOrbitDetermination.jl")
println("    ✓ Exported official CCSDS 508.0-B-1 CDM: $cdm_file")
println("\n================================================================================")
println(" DEMONSTRATION COMPLETE: READY FOR ISRO / FLIGHT DYNAMICS REVIEW")
println("================================================================================")
