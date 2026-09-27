# Orbit Propagator & Differential Equations Integrator Interface

struct OrbitPropagatorOptions
    include_j2::Bool
    include_drag::Bool
    include_srp::Bool
    Cd::Float64
    Cr::Float64
    area_to_mass::Float64 # m^2/kg
    reltol::Float64
    abstol::Float64
end

function OrbitPropagatorOptions(;
    include_j2::Bool=true,
    include_drag::Bool=false,
    include_srp::Bool=false,
    Cd::Float64=2.2,
    Cr::Float64=1.2,
    area_to_mass::Float64=0.01,
    reltol::Float64=1e-8,
    abstol::Float64=1e-8
)
    return OrbitPropagatorOptions(include_j2, include_drag, include_srp, Cd, Cr, area_to_mass, reltol, abstol)
end

"""
    orbit_dynamics(u, p, t)

ODE state derivative for orbit propagation: state u = [x, y, z, vx, vy, vz].
Parameters p = (opts::OrbitPropagatorOptions, mu::Float64)
Guarantees concrete element type return matching input vector u.
"""
function orbit_dynamics(u, p, t)
    opts, mu = p
    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    # Central body acceleration
    a_tot = two_body_acceleration(r, mu)

    # Add perturbations as requested
    if opts.include_j2
        a_tot += j2_perturbation(r; mu=mu)
    end
    if opts.include_drag
        a_tot += drag_perturbation(r, v, opts.Cd, opts.area_to_mass)
    end
    if opts.include_srp
        a_tot += srp_perturbation(r, opts.Cr, opts.area_to_mass)
    end

    T = eltype(u)
    return T[v[1], v[2], v[3], a_tot[1], a_tot[2], a_tot[3]]
end

"""
    propagate_orbit(u0::AbstractVector, tspan::Tuple{Float64, Float64}, saveat::Vector{Float64}, opts::OrbitPropagatorOptions=OrbitPropagatorOptions(); mu::Float64=MU_EARTH)

Propagate Cartesian orbit state u0 over tspan and save solution at saveat times.
Guarantees concrete floating point state type promotion for Tsit5 solver compatibility.
"""
function propagate_orbit(u0::AbstractVector, tspan::Tuple{Float64, Float64}, saveat::Vector{Float64}, opts::OrbitPropagatorOptions=OrbitPropagatorOptions(); mu::Float64=MU_EARTH)
    T = typeof(float(u0[1]))
    u0_concrete = T[float(x) for x in u0]
    prob = ODEProblem(orbit_dynamics, u0_concrete, tspan, (opts, mu))
    sol = solve(prob, Tsit5(), saveat=saveat, reltol=opts.reltol, abstol=opts.abstol)
    return sol
end
