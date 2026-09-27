# High-Order Spherical Harmonics Earth Gravity Model (J2, J3, J4, J5 Zonal Harmonics)

const J3_EARTH = -2.5324105e-6
const J4_EARTH = -1.6199216e-6
const J5_EARTH = -2.2788820e-7

"""
    zonal_harmonic_perturbation(r::SVector{3, T}; mu::Float64=MU_EARTH, R_e::Float64=R_EARTH) where T

Compute acceleration vector due to high-order Earth zonal harmonics (J2, J3, J4, J5) in ECI frame.
"""
function zonal_harmonic_perturbation(r::SVector{3, T}; mu::Float64=MU_EARTH, R_e::Float64=R_EARTH) where T
    x, y, z = r[1], r[2], r[3]
    r_norm = norm(r)
    z_r = z / r_norm

    # J2 acceleration
    a_j2 = j2_perturbation(r; mu=mu, R_e=R_e, J2=J2_EARTH)

    # J3 acceleration
    f_j3 = 2.5 * J3_EARTH * mu * (R_e^3) / (r_norm^6)
    ax_j3 = f_j3 * x * (7.0 * z_r^3 - 3.0 * z_r)
    ay_j3 = f_j3 * y * (7.0 * z_r^3 - 3.0 * z_r)
    az_j3 = f_j3 * (3.5 * z_r^4 - 3.0 * z_r^2 + 0.5) * r_norm
    a_j3 = SVector{3, T}(ax_j3, ay_j3, az_j3)

    # J4 acceleration
    f_j4 = 0.625 * J4_EARTH * mu * (R_e^4) / (r_norm^7)
    ax_j4 = f_j4 * x * (35.0 * z_r^4 - 30.0 * z_r^2 + 3.0)
    ay_j4 = f_j4 * y * (35.0 * z_r^4 - 30.0 * z_r^2 + 3.0)
    az_j4 = f_j4 * z * (35.0 * z_r^4 - 50.0 * z_r^2 + 15.0)
    a_j4 = SVector{3, T}(ax_j4, ay_j4, az_j4)

    return a_j2 + a_j3 + a_j4
end
