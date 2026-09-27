# Orbital Perturbations: J2 Oblateness, Atmospheric Drag, and Solar Radiation Pressure (SRP)

"""
    j2_perturbation(r::SVector{3, T}; mu::Float64=MU_EARTH, R_e::Float64=R_EARTH, J2::Float64=J2_EARTH) where T

Compute acceleration due to Earth J2 oblateness effect in ECI frame.
"""
function j2_perturbation(r::SVector{3, T}; mu::Float64=MU_EARTH, R_e::Float64=R_EARTH, J2::Float64=J2_EARTH) where T
    x, y, z = r[1], r[2], r[3]
    r_norm = norm(r)
    r_eff = r_norm < 6000.0 ? T(6000.0) : r_norm
    factor = 1.5 * J2 * mu * (R_e^2) / (r_eff^5)
    z_sq_ratio = (z / r_eff)^2

    ax = factor * x * (5.0 * z_sq_ratio - 1.0)
    ay = factor * y * (5.0 * z_sq_ratio - 1.0)
    az = factor * z * (5.0 * z_sq_ratio - 3.0)

    return SVector{3, T}(ax, ay, az)
end

"""
    atmospheric_density_exponential(alt_km::T) where T

Simple exponential atmospheric density model.
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
    srp_perturbation(r::SVector{3, T}, Cr::Float64, area_to_mass::Float64) where T

Compute Solar Radiation Pressure acceleration (km/s^2).
"""
function srp_perturbation(r::SVector{3, T}, Cr::Float64, area_to_mass::Float64) where T
    s_sun = SVector{3, T}(1.0, 0.0, 0.0)
    P_sun_km = 4.56e-6 * 1e-3
    A_m_km2 = area_to_mass * 1e-6
    
    proj = dot(r, s_sun)
    if proj < 0.0 && norm(r - proj * s_sun) < R_EARTH
        return SVector{3, T}(0.0, 0.0, 0.0)
    end
    
    return Cr * P_sun_km * A_m_km2 * s_sun
end
