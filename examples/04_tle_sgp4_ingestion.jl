# Tutorial 04: Real TLE Parsing, SGP4 Analytical Propagation, and Bayesian OD Refinement

using BayesianOrbitDetermination
using Printf

println("=== BayesianOrbitDetermination.jl Tutorial 04 ===")
println("Target: ISS TLE Parsing & SGP4 Dynamics Ingestion\n")

tle_line1 = "1 25544U 98067A   26270.50000000  .00016717  00000-0  30000-3 0  9993"
tle_line2 = "2 25544  51.6400 208.1234 0004500  65.1234 295.0000 15.49500000420002"

tle = parse_tle(tle_line1, tle_line2, "ISS (ZARYA)")
@printf("Parsed Satellite: %s (NORAD ID: %d)\n", tle.sat_name, tle.norad_id)
@printf("Inclination: %.2f deg | Eccentricity: %.5f | RAAN: %.2f deg\n", tle.inc_deg, tle.e, tle.raan_deg)

c_sgp4 = propagate_sgp4(tle, 1800.0) # 30 min propagation
@printf("SGP4 Position at +30 min (km): [%.2f, %.2f, %.2f]\n", c_sgp4.r[1], c_sgp4.r[2], c_sgp4.r[3])

println("\nTutorial 04 complete successfully!")
