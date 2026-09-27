# SGP4 / TLE Analytical Orbital Element Propagator & Parser

struct TLE
    sat_name::String
    norad_id::Int
    epoch_year::Int
    epoch_day::Float64
    bstar::Float64
    inc_deg::Float64
    raan_deg::Float64
    e::Float64
    argp_deg::Float64
    mean_anomaly_deg::Float64
    mean_motion_rev_day::Float64
end

"""
    parse_tle(line1::String, line2::String, name::String="SATELLITE") -> TLE

Parse standard TLE Line 1 and Line 2 into structured `TLE` type.
"""
function parse_tle(line1::String, line2::String, name::String="SATELLITE")
    norad_id = parse(Int, strip(line1[3:7]))
    epoch_yr = parse(Int, line1[19:20])
    epoch_day = parse(Float64, line1[21:32])
    
    bstar_str = strip(line1[54:61])
    bstar_val = parse(Float64, bstar_str[1:end-2]) * 10.0^(parse(Float64, bstar_str[end-1:end])) * 1e-5

    inc = parse(Float64, strip(line2[9:16]))
    raan = parse(Float64, strip(line2[18:25]))
    e = parse(Float64, "0." * strip(line2[27:33]))
    argp = parse(Float64, strip(line2[35:42]))
    ma = parse(Float64, strip(line2[44:51]))
    mm = parse(Float64, strip(line2[53:63]))

    return TLE(name, norad_id, epoch_yr, epoch_day, bstar_val, inc, raan, e, argp, ma, mm)
end

"""
    tle_to_keplerian(tle::TLE) -> KeplerianState

Convert TLE parameters to instantaneous osculating KeplerianState.
"""
function tle_to_keplerian(tle::TLE)
    # Mean motion (rev/day) to semi-major axis a (km)
    n_rad_s = tle.mean_motion_rev_day * (2pi / 86400.0)
    a = (MU_EARTH / (n_rad_s^2))^(1/3)
    
    return KeplerianState(a, tle.e, deg2rad(tle.inc_deg), deg2rad(tle.raan_deg), deg2rad(tle.argp_deg), deg2rad(tle.mean_anomaly_deg))
end

"""
    propagate_sgp4(tle::TLE, t_sec::Float64) -> CartesianState

Analytical SGP4 propagation approximation over time delta `t_sec` from TLE epoch.
"""
function propagate_sgp4(tle::TLE, t_sec::Float64)
    k_init = tle_to_keplerian(tle)
    n_rad_s = tle.mean_motion_rev_day * (2pi / 86400.0)
    
    # Propagate mean anomaly and secular J2 RAAN drift
    delta_M = n_rad_s * t_sec
    nu_new = mod(k_init.nu + delta_M, 2pi)

    # Secular J2 drift rate for RAAN (rad/s)
    p = k_init.a * (1.0 - k_init.e^2)
    raan_dot = -1.5 * J2_EARTH * (R_EARTH / p)^2 * n_rad_s * cos(k_init.i)
    raan_new = mod(k_init.RAAN + raan_dot * t_sec, 2pi)

    k_t = KeplerianState(k_init.a, k_init.e, k_init.i, raan_new, k_init.argp, nu_new)
    return kepler_to_cartesian(k_t)
end
