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
    include_j2::Bool = true,
    include_drag::Bool = false,
    include_srp::Bool = false,
    Cd::Float64 = 2.2,
    Cr::Float64 = 1.2,
    area_to_mass::Float64 = 0.01,
    reltol::Float64 = 1e-12,
    abstol::Float64 = 1e-12,
)
    return OrbitPropagatorOptions(
        include_j2, include_drag, include_srp, Cd, Cr, area_to_mass, reltol, abstol,
    )
end

"""
    orbit_dynamics(u, p, t)

ODE state derivative for orbit propagation.

State: `u = [x, y, z, vx, vy, vz]`
Parameters: `p = (opts::OrbitPropagatorOptions, mu::Float64)`

Returns a vector with the same element type as `u`, so this works with
ForwardDiff / Turing dual numbers.
"""
function orbit_dynamics(u, p, t)
    opts, mu = p
    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    a_tot = two_body_acceleration(r, mu)

    if opts.include_j2
        a_tot += j2_perturbation(r; mu = mu)
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
    propagate_orbit(u0, tspan, saveat, opts = OrbitPropagatorOptions();
                    mu = MU_EARTH)

Propagate a Cartesian orbit state `u0` over `tspan` and save at `saveat` times.

The element type of the returned solution matches the element type of `u0`,
so this function can be called with `ForwardDiff.Dual` numbers inside a
Turing model.
"""
function propagate_orbit(
    u0::AbstractVector,
    tspan::Tuple{Float64, Float64},
    saveat::Vector{Float64},
    opts::OrbitPropagatorOptions = OrbitPropagatorOptions();
    mu::Float64 = MU_EARTH,
)
    u0_concrete = collect(u0)
    prob = ODEProblem(orbit_dynamics, u0_concrete, tspan, (opts, mu))
    sol = solve(
        prob,
        Tsit5();
        saveat = saveat,
        reltol = opts.reltol,
        abstol = opts.abstol,
    )
    return sol
end