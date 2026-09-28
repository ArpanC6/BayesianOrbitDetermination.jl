# Two-body Keplerian Dynamics and Coordinate Transformations

struct KeplerianState
    a::Float64 # Semi-major axis (km)
    e::Float64 # Eccentricity
    i::Float64 # Inclination (rad)
    RAAN::Float64 # Right Ascension of Ascending Node (rad)
    argp::Float64 # Argument of Perigee (rad)
    nu::Float64 # True Anomaly (rad)
end

struct CartesianState{T<:Real}
    r::SVector{3, T} # Position vector [x, y, z] (km)
    v::SVector{3, T} # Velocity vector [vx, vy, vz] (km/s)
end

"""
    two_body_acceleration(r::SVector{3, T}, mu::Float64=MU_EARTH) where T

Compute standard central-body gravitational acceleration with smooth singularity
regularization to prevent ODE integrator instability during MCMC parameter exploration.
The soft-floor ensures r_eff is never below ~R_EARTH while remaining differentiable
everywhere (critical for ForwardDiff/NUTS HMC).
"""
function two_body_acceleration(r::SVector{3, T}, mu::Float64=MU_EARTH) where T
    r_norm = norm(r)
    r_eff = sqrt(r_norm^2 + R_EARTH^2 * exp(-r_norm / R_EARTH))
    return -mu * r / (r_eff^3)
end

"""
    kepler_to_cartesian(k::KeplerianState; mu::Float64=MU_EARTH) -> CartesianState{Float64}

Convert Keplerian orbital elements to Cartesian inertial state (ECI).
"""
function kepler_to_cartesian(k::KeplerianState; mu::Float64=MU_EARTH)
    p = k.a * (1.0 - k.e^2)
    r_pqw = SVector(p * cos(k.nu) / (1.0 + k.e * cos(k.nu)),
                    p * sin(k.nu) / (1.0 + k.e * cos(k.nu)),
                    0.0)
    v_pqw = SVector(-sqrt(mu / p) * sin(k.nu),
                     sqrt(mu / p) * (k.e + cos(k.nu)),
                     0.0)

    # Rotation matrix PQW -> ECI
    c_node, s_node = cos(k.RAAN), sin(k.RAAN)
    c_inc, s_inc = cos(k.i), sin(k.i)
    c_arg, s_arg = cos(k.argp), sin(k.argp)

    R = SMatrix{3, 3}([
        (c_node * c_arg - s_node * s_arg * c_inc)  (-c_node * s_arg - s_node * c_arg * c_inc)  (s_node * s_inc);
        (s_node * c_arg + c_node * s_arg * c_inc)  (-s_node * s_arg + c_node * c_arg * c_inc)  (-c_node * s_inc);
        (s_arg * s_inc)                            (c_arg * s_inc)                             (c_inc)
    ])

    return CartesianState(R * r_pqw, R * v_pqw)
end

"""
    cartesian_to_kepler(c::CartesianState; mu::Float64=MU_EARTH) -> KeplerianState

Convert Cartesian inertial state (ECI) to Keplerian orbital elements.
"""
function cartesian_to_kepler(c::CartesianState; mu::Float64=MU_EARTH)
    r_vec = c.r
    v_vec = c.v
    r = norm(r_vec)
    v = norm(v_vec)

    h_vec = cross(r_vec, v_vec)
    h = norm(h_vec)

    n_vec = cross(SVector(0.0, 0.0, 1.0), h_vec)
    n = norm(n_vec)

    e_vec = ((v^2 - mu / r) * r_vec - dot(r_vec, v_vec) * v_vec) / mu
    e = norm(e_vec)

    energy = 0.5 * v^2 - mu / r
    a = -mu / (2.0 * energy)

    inc = acos(clamp(h_vec[3] / h, -1.0, 1.0))

    RAAN = n != 0 ? acos(clamp(n_vec[1] / n, -1.0, 1.0)) : 0.0
    if n_vec[2] < 0
        RAAN = 2pi - RAAN
    end

    argp = (n != 0 && e > 1e-11) ? acos(clamp(dot(n_vec, e_vec) / (n * e), -1.0, 1.0)) : 0.0
    if e_vec[3] < 0
        argp = 2pi - argp
    end

    nu = e > 1e-11 ? acos(clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0)) : 0.0
    if dot(r_vec, v_vec) < 0
        nu = 2pi - nu
    end

    return KeplerianState(a, e, inc, RAAN, argp, nu)
end