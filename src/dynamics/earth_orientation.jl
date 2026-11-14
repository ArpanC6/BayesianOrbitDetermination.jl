# Earth orientation helpers for real tracking data.

"""
    gmst_rad(jd_utc) -> Float64

Greenwich mean sidereal time in radians from the IAU 1982 expression. UTC is
used in place of UT1; the difference is below one second of time, which is
negligible for the range geometry of a satellite in low Earth orbit.
"""
function gmst_rad(jd_utc::Float64)
    d = jd_utc - 2451545.0
    T = d / 36525.0
    theta = 280.46061837 + 360.98564736629 * d + 0.000387933 * T^2 - T^3 / 38710000.0
    return deg2rad(mod(theta, 360.0))
end

"""
    itrf_to_eci_state(r_itrf, v_itrf, jd_utc; omega_earth = OMEGA_EARTH)

Rotate an ITRF position and velocity (km, km/s) into the inertial frame of
date. The velocity includes the Earth-rotation term `omega x r`.
Returns `(r_eci, v_eci, theta)` with `theta = gmst_rad(jd_utc)`.
"""
function itrf_to_eci_state(r_itrf::SVector{3, Float64}, v_itrf::SVector{3, Float64},
                           jd_utc::Float64; omega_earth::Float64 = OMEGA_EARTH)
    theta = gmst_rad(jd_utc)
    c, s = cos(theta), sin(theta)
    R = SMatrix{3, 3, Float64}([c -s 0.0; s c 0.0; 0.0 0.0 1.0])
    omega = SVector{3, Float64}(0.0, 0.0, omega_earth)
    r_eci = R * r_itrf
    v_eci = R * (v_itrf + cross(omega, r_itrf))
    return r_eci, v_eci, theta
end
