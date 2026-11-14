# Orbital Perturbations: J2 Oblateness, Atmospheric Drag, and Solar Radiation Pressure (SRP)

"""
    j2_perturbation(r::SVector{3, T}; mu::Float64=MU_EARTH, R_e::Float64=R_EARTH, J2::Float64=J2_EARTH) where T

Compute acceleration due to Earth J2 oblateness effect in ECI frame.
"""
function j2_perturbation(r::SVector{3, T}; mu::Float64=MU_EARTH, R_e::Float64=R_EARTH, J2::Float64=J2_EARTH) where T
    x, y, z = r[1], r[2], r[3]
    r_norm = norm(r)
    factor = 1.5 * J2 * mu * (R_e^2) / (r_norm^5)
    z_sq_ratio = (z / r_norm)^2

    ax = factor * x * (5.0 * z_sq_ratio - 1.0)
    ay = factor * y * (5.0 * z_sq_ratio - 1.0)
    az = factor * z * (5.0 * z_sq_ratio - 3.0)

    return SVector{3, T}(ax, ay, az)
end

"""
    atmospheric_density_exponential(alt_km::T) where T

Simple exponential atmospheric density model. Returns density in kg/km^3,
which is the unit that combines with `area_to_mass` in km^2/kg in
[`drag_perturbation`](@ref) to give acceleration in km/s^2.
"""
function atmospheric_density_exponential(alt_km::T) where T
    if alt_km < 0.0
        return T(1.225e9)
    end
    rho0 = 2.789e-10 * 1e9
    h0 = 200.0
    H = 50.0
    return rho0 * exp(-(alt_km - h0) / H)
end

"""
    drag_perturbation(r::SVector{3, T}, v::SVector{3, T}, Cd::Float64, area_to_mass::Float64; omega_earth::Float64=OMEGA_EARTH) where T

Compute atmospheric drag acceleration (km/s^2).
"""
function drag_perturbation(r::SVector{3, T}, v::SVector{3, T}, Cd::Float64, area_to_mass::Float64; omega_earth::Float64=OMEGA_EARTH) where T
    alt = norm(r) - R_EARTH
    if alt > 1000.0
        return SVector{3, T}(0.0, 0.0, 0.0)
    end
    rho = atmospheric_density_exponential(alt)
    
    v_rel = v - cross(SVector{3, T}(0.0, 0.0, omega_earth), r)
    v_rel_norm = norm(v_rel)
    A_m_km2 = area_to_mass * 1e-6
    
    a_drag = -0.5 * Cd * A_m_km2 * rho * v_rel_norm * v_rel
    return a_drag
end

"""
    sun_direction_eci(jd_utc) -> SVector{3, Float64}

Unit vector from the Earth toward the Sun in the mean equator and equinox of
date, from the low-precision solar ephemeris in Vallado (2013). Accurate to
about 0.01 deg, which is ample for solar radiation pressure.
"""
function sun_direction_eci(jd_utc::Float64)
    T = (jd_utc - 2451545.0) / 36525.0
    L = deg2rad(mod(280.460 + 36000.771 * T, 360.0))
    M = deg2rad(mod(357.5291092 + 35999.05034 * T, 360.0))
    lambda = L + deg2rad(1.914666471) * sin(M) + deg2rad(0.019994643) * sin(2M)
    eps = deg2rad(23.439291 - 0.0130042 * T)
    return SVector(cos(lambda), cos(eps) * sin(lambda), sin(eps) * sin(lambda))
end

"""
    sun_distance_au(jd_utc) -> Float64

Earth-Sun distance in astronomical units from the same low-precision
ephemeris as [`sun_direction_eci`](@ref).
"""
function sun_distance_au(jd_utc::Float64)
    T = (jd_utc - 2451545.0) / 36525.0
    M = deg2rad(mod(357.5291092 + 35999.05034 * T, 360.0))
    return 1.000140612 - 0.016708617 * cos(M) - 0.000139589 * cos(2M)
end

"""
    srp_perturbation(r, t, Cr, area_to_mass; epoch_jd = 2451545.0)

Compute Solar Radiation Pressure acceleration (km/s^2) at `t` seconds from the
epoch whose Julian date is `epoch_jd`. The Sun direction follows the
low-precision ephemeris and the acceleration is scaled by the inverse square
of the Earth-Sun distance. A cylindrical shadow model zeroes the acceleration
when the satellite is behind the Earth.

`area_to_mass` is in m^2/kg. The solar pressure at 1 AU is 4.56e-6 N/m^2.
"""
function srp_perturbation(r::SVector{3, T}, t, Cr::Float64,
                          area_to_mass::Float64; epoch_jd::Float64 = 2451545.0) where T
    jd = epoch_jd + t / 86400.0
    s_sun = SVector{3, T}(sun_direction_eci(jd))
    au_ratio = sun_distance_au(jd)
    P_sun = 4.56e-6 / au_ratio^2          # N/m^2
    a_srp = Cr * P_sun * area_to_mass * 1e-3  # km/s^2

    proj = dot(r, s_sun)
    if proj < 0.0 && norm(r - proj * s_sun) < R_EARTH
        return SVector{3, T}(0.0, 0.0, 0.0)
    end

    return a_srp * s_sun
end
