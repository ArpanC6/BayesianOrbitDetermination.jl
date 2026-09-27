# Space Traffic Management: Conjunction Assessment & 2D Collision Probability (Pc) Engine

struct ConjunctionEvent
    sat1_name::String
    sat2_name::String
    t_tca_sec::Float64        # Time of Closest Approach
    r1_eci::SVector{3, Float64}
    v1_eci::SVector{3, Float64}
    cov1_eci::Matrix{Float64}
    r2_eci::SVector{3, Float64}
    v2_eci::SVector{3, Float64}
    cov2_eci::Matrix{Float64}
    hard_body_radius_m::Float64 # Hard body sphere combined radius (m)
end

"""
    compute_collision_probability_foster(conj::ConjunctionEvent) -> Float64

Compute 2D Foster-Elrod Probability of Collision (Pc) in the encounter plane.
Project 3D position covariances into the 2D plane perpendicular to relative velocity vector v_rel.
"""
function compute_collision_probability_foster(conj::ConjunctionEvent)
    r_rel = conj.r2_eci - conj.r1_eci # km
    v_rel = conj.v2_eci - conj.v1_eci # km/s
    v_rel_norm = norm(v_rel)

    if v_rel_norm < 1e-6
        return 0.0
    end

    # Encounter frame axes: y_enc along v_rel, z_enc along r_rel x v_rel, x_enc completes right-hand frame
    y_enc = v_rel / v_rel_norm
    z_enc = cross(r_rel, v_rel)
    z_enc = norm(z_enc) > 1e-9 ? z_enc / norm(z_enc) : SVector(0.0, 0.0, 1.0)
    x_enc = cross(y_enc, z_enc)

    R_enc = SMatrix{3, 3}([x_enc'; y_enc'; z_enc']) # Transformation matrix ECI -> Encounter

    # Combined 3D covariance matrix in ECI
    C_combined_eci = conj.cov1_eci[1:3, 1:3] + conj.cov2_eci[1:3, 1:3]
    C_enc = R_enc * C_combined_eci * R_enc'

    # Project into 2D Encounter Plane (x_enc, z_enc)
    C_2d = [C_enc[1, 1] C_enc[1, 3]; C_enc[3, 1] C_enc[3, 3]]
    r_2d = [dot(r_rel, x_enc), dot(r_rel, z_enc)] # Relative miss vector in 2D plane

    # Hard body radius km
    r_hb_km = conj.hard_body_radius_m * 1e-3

    inv_C2d = try
        inv(C_2d)
    catch
        pinv(C_2d)
    end
    det_C2d = det(C_2d)

    if det_C2d <= 0.0
        return 0.0
    end

    # Numerical polar integration over combined hard body disk radius r_hb
    n_r, n_theta = 40, 72
    dr = r_hb_km / n_r
    dtheta = 2pi / n_theta
    pc_sum = 0.0

    factor = 1.0 / (2pi * sqrt(det_C2d))

    for i in 1:n_r
        r_val = (i - 0.5) * dr
        for j in 1:n_theta
            theta_val = (j - 0.5) * dtheta
            dx = r_2d[1] + r_val * cos(theta_val)
            dz = r_2d[2] + r_val * sin(theta_val)
            vec = [dx, dz]

            mahalanobis_sq = dot(vec, inv_C2d * vec)
            pdf_val = factor * exp(-0.5 * mahalanobis_sq)
            pc_sum += pdf_val * r_val * dr * dtheta
        end
    end

    return pc_sum
end
