# TLE parsing and SGP4 propagation.
#
# The numerical propagator in this package is not SGP4. SGP4 is a specific
# analytical theory for TLEs, and it is delegated to SatelliteToolboxSgp4.jl
# rather than approximated with Keplerian motion plus a J2 drift.

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
    line1::String
    line2::String
end

# Kept for callers that built a TLE from elements; SGP4 needs the raw lines.
TLE(name, norad_id, epoch_year, epoch_day, bstar, inc, raan, e, argp, ma, mm) =
    TLE(name, norad_id, epoch_year, epoch_day, bstar, inc, raan, e, argp, ma, mm, "", "")

"""
    parse_tle(line1::String, line2::String, name::String="SATELLITE") -> TLE

Parse standard TLE Line 1 and Line 2 into a structured `TLE`. The raw lines are
kept so that [`propagate_sgp4`](@ref) can delegate to a real SGP4 implementation.
"""
function parse_tle(line1::String, line2::String, name::String = "SATELLITE")
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

    return TLE(name, norad_id, epoch_yr, epoch_day, bstar_val, inc, raan, e, argp,
               ma, mm, line1, line2)
end

"""
    tle_to_keplerian(tle::TLE) -> KeplerianState

Osculating Keplerian elements at the TLE epoch, obtained by evaluating SGP4 at
`t = 0` and converting the resulting TEME state. This is not the same as the
mean elements stored in the TLE.
"""
function tle_to_keplerian(tle::TLE)
    return cartesian_to_kepler(propagate_sgp4(tle, 0.0))
end

"""
    propagate_sgp4(tle::TLE, t_sec::Float64; verify_checksum::Bool=false) -> CartesianState

Propagate a TLE with the real SGP4/SDP4 model from `SatelliteToolboxSgp4.jl`.

`t_sec` is seconds from the TLE epoch. The returned position (km) and velocity
(km/s) are in the TEME frame, which is what SGP4 is defined in; convert TEME to
EME2000 before mixing these states with the package's inertial propagator.

TLE checksums are not verified by default because example TLEs in tutorials and
tests often use placeholder checksums. Pass `verify_checksum = true` for real
distributed TLEs.
"""
function propagate_sgp4(tle::TLE, t_sec::Float64; verify_checksum::Bool = false)
    (isempty(tle.line1) || isempty(tle.line2)) &&
        throw(ArgumentError("TLE has no raw lines; SGP4 needs line1 and line2. Use parse_tle."))

    parsed = SatelliteToolboxTle.read_tle(tle.line1, tle.line2;
                                          name = tle.sat_name,
                                          verify_checksum = verify_checksum)
    sgp4d = SatelliteToolboxSgp4.sgp4_init(parsed)
    r, v = SatelliteToolboxSgp4.sgp4!(sgp4d, t_sec / 60.0)
    return CartesianState(SVector{3}(r...), SVector{3}(v...))
end
