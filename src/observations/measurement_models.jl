# Ground Station Measurement Models: Range, Range Rate (Doppler), Azimuth, Elevation
#
# This module provides:
# - Geodetic to ECEF conversion (WGS-84 ellipsoid)
# - ECEF to ECI rotation (Earth rotation)
# - Satellite-to-station geometry (range, Doppler, azimuth, elevation)
# - Synthetic tracking data generation with visibility filtering
#
# All functions are generic over the number type (Float64, ForwardDiff.Dual)
# so they can be used inside Turing models without modification.

const C_LIGHT_KM_S = 299792.458  # speed of light (km/s)

"""
    GroundStation

Geodetic coordinates of a ground-based tracking station.

# Fields
- `name::String`       : Station identifier
- `lat_deg::Float64`   : Geodetic latitude (degrees)
- `lon_deg::Float64`   : Geodetic longitude (degrees)
- `alt_km::Float64`    : Altitude above WGS-84 ellipsoid (km)
"""
struct GroundStation
    name::String
    lat_deg::Float64
    lon_deg::Float64
    alt_km::Float64
end

"""
    ground_station_ecef(gs::GroundStation) -> SVector{3, Float64}

Convert geodetic coordinates (lat, lon, alt) to ECEF Cartesian (km).

Uses the WGS-84 ellipsoid with eccentricity squared e² = 0.00669437999014.
"""
function ground_station_ecef(gs::GroundStation)
    phi = deg2rad(gs.lat_deg)
    lambda = deg2rad(gs.lon_deg)
    h = gs.alt_km

    e2 = 0.00669437999014
    N = R_EARTH / sqrt(1.0 - e2 * sin(phi)^2)

    x = (N + h) * cos(phi) * cos(lambda)
    y = (N + h) * cos(phi) * sin(lambda)
    z = (N * (1.0 - e2) + h) * sin(phi)

    return SVector(x, y, z)
end

"""
    ecef_to_eci(r_ecef::SVector{3, Float64}, t_sec::Float64;
                omega_earth::Float64 = OMEGA_EARTH) -> SVector{3, Float64}

Rotate an ECEF vector to ECI frame at time `t_sec` from epoch.

# Arguments
- `r_ecef`     : Position vector in ECEF (km)
- `t_sec`      : Time from epoch (seconds)
- `omega_earth`: Earth rotation rate (rad/s), default `OMEGA_EARTH`

# Returns
Position vector in ECI (km).
"""
function ecef_to_eci(r_ecef::SVector{3, Float64}, t_sec::Float64;
                     omega_earth::Float64 = OMEGA_EARTH)
    theta = omega_earth * t_sec
    c, s = cos(theta), sin(theta)
    R_z = SMatrix{3, 3}([c -s 0.0; s c 0.0; 0.0 0.0 1.0])
    return R_z * r_ecef
end

"""
    eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t_sec;
                                       omega_earth = OMEGA_EARTH,
                                       light_time = true)

Compute topocentric Azimuth, Elevation, Range, and Range Rate (Doppler)
from a satellite state and a ground station.

# Arguments
- `r_sat`  : Satellite ECI position (km)
- `v_sat`  : Satellite ECI velocity (km/s)
- `gs`     : `GroundStation`
- `t_sec`  : Time from epoch (seconds)
- `light_time`: solve the two-way light-time equation (default). Set to
  `false` for the instantaneous geometric range and rate, e.g. when
  comparing against a model that quotes geometry at the reception epoch.

# Returns
Named tuple `(azimuth, elevation, range, range_rate)`.

# Light time
A range measurement is the two-way time of flight, so the satellite is at
the bounce epoch `t_sec - tau_down`, not at the reception epoch. At a
1300 km range the one-way light time is 4 ms, which is tens of metres of
along-track motion for a LEO satellite. The down-leg and up-leg light
times are solved by fixed-point iteration, with the satellite displaced by
`v_sat * tau` and the station rotated back through `omega_earth * tau`.
This is the dominant correction for SLR data; ignoring it biases the fit
by ~20 m.

This function is generic over the number type `T` (Float64, Dual, etc.),
so it can be used inside automatic-differentiation pipelines.
"""
function eci_to_station_azel_range_doppler(
    r_sat::SVector{3, T},
    v_sat::SVector{3, T},
    gs::GroundStation,
    t_sec::Float64;
    omega_earth::Float64 = OMEGA_EARTH,
    light_time::Bool = true,
) where {T}
    r_gs_ecef = ground_station_ecef(gs)
    theta = omega_earth * t_sec
    ct, st = cos(theta), sin(theta)
    xe, ye, ze = r_gs_ecef
    r_gs_eci = SVector{3, T}(ct * xe - st * ye, st * xe + ct * ye, ze)
    v_gs_eci = SVector{3, T}(-omega_earth * r_gs_eci[2],
                             omega_earth * r_gs_eci[1], zero(T))

    tau_down = zero(T)
    tau_up = zero(T)
    if light_time
        for _ in 1:3
            sat_b = r_sat - v_sat * tau_down
            rho_down = norm(sat_b - r_gs_eci)
            tau_down = rho_down / C_LIGHT_KM_S

            theta_up = theta - omega_earth * (tau_down + tau_up)
            cu, su = cos(theta_up), sin(theta_up)
            r_gs_up = SVector{3, T}(cu * xe - su * ye, su * xe + cu * ye, ze)
            tau_up = norm(sat_b - r_gs_up) / C_LIGHT_KM_S
        end
        # The observable is half the round-trip time of flight.
        range_km = C_LIGHT_KM_S * (tau_up + tau_down) / 2
        theta_b = theta - omega_earth * tau_down
        cb, sb = cos(theta_b), sin(theta_b)
        r_gs_b = SVector{3, T}(cb * xe - sb * ye, sb * xe + cb * ye, ze)
        v_gs_b = SVector{3, T}(-omega_earth * r_gs_b[2],
                               omega_earth * r_gs_b[1], zero(T))
        rho_vec = (r_sat - v_sat * tau_down) - r_gs_b
        rho_dot_vec = v_sat - v_gs_b
        lambda = deg2rad(gs.lon_deg) + theta_b
    else
        rho_vec = r_sat - r_gs_eci
        rho_dot_vec = v_sat - v_gs_eci
        range_km = norm(rho_vec)
        lambda = deg2rad(gs.lon_deg) + theta
    end

    range_rate = dot(rho_vec, rho_dot_vec) / norm(rho_vec)

    phi = deg2rad(gs.lat_deg)
    sin_lat, cos_lat = sin(phi), cos(phi)
    sin_lon, cos_lon = sin(lambda), cos(lambda)

    R_sez = SMatrix{3, 3, T}([
        (sin_lat * cos_lon) (sin_lat * sin_lon) (-cos_lat)
        (-sin_lon)          (cos_lon)           (zero(T))
        (cos_lat * cos_lon) (cos_lat * sin_lon) (sin_lat)
    ])

    rho_sez = R_sez * rho_vec
    S, E, Z = rho_sez[1], rho_sez[2], rho_sez[3]

    elevation = asin(clamp(Z / norm(rho_vec), -1.0, 1.0))
    azimuth = atan(E, -S)
    if azimuth < 0.0
        azimuth += 2pi
    end

    return (azimuth = azimuth, elevation = elevation,
            range = range_km, range_rate = range_rate)
end

"""
    simulate_tracking_data(u0, t_obs, gs, sigma_range, sigma_doppler;
                            opts = OrbitPropagatorOptions(),
                            min_elevation_deg = 5.0)

Generate synthetic ground-station tracking data (range and Doppler)
for a satellite over the requested observation times.

# Arguments
- `u0`               : Initial satellite state `[x, y, z, vx, vy, vz]` (km, km/s) at epoch `t = 0`
- `t_obs`            : Observation times (seconds from epoch 0)
- `gs`               : Ground station
- `sigma_range`      : Range measurement noise standard deviation (km)
- `sigma_doppler`    : Doppler measurement noise standard deviation (km/s)
- `opts`             : Propagation options
- `min_elevation_deg`: Minimum elevation for a valid pass (degrees)

# Returns
Named tuple `(t_obs, range, doppler, sol)`.

Observations are filtered to those where the satellite is above the
horizon (`elevation > min_elevation_deg`). If none pass, empty vectors are
returned; callers must treat that as a failed scenario rather than fitting
below-horizon data.
"""
function simulate_tracking_data(
    u0::Vector{Float64},
    t_obs::Vector{Float64},
    gs::GroundStation,
    sigma_range::Float64,
    sigma_doppler::Float64;
    opts = OrbitPropagatorOptions(),
    min_elevation_deg::Float64 = 5.0,
)
    tspan = (0.0, t_obs[end])
    sol = propagate_orbit(u0, tspan, t_obs, opts)

    min_elev = deg2rad(min_elevation_deg)

    valid_t = Float64[]
    obs_range = Float64[]
    obs_doppler = Float64[]

    for (i, t) in enumerate(t_obs)
        r_sat = SVector{3, Float64}(sol.u[i][1:3])
        v_sat = SVector{3, Float64}(sol.u[i][4:6])
        meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t)

        if meas.elevation > min_elev
            push!(valid_t, t)
            push!(obs_range, meas.range + randn() * sigma_range)
            push!(obs_doppler, meas.range_rate + randn() * sigma_doppler)
        end
    end

    return (t_obs = valid_t, range = obs_range,
            doppler = obs_doppler, sol = sol)
end

"""
    select_pass_times(u0, gs, opts; t_end = 24 * 3600.0, dt = 10.0,
                      min_elevation_deg = 10.0, per_pass = 8, max_passes = 3)

Return observation times for the first `max_passes` visibility windows of a
satellite above `min_elevation_deg`, with `per_pass` samples spread over each
window. Sparse tracking means a few short passes, not a continuous arc.

Returns `(times, n_passes)`; `times` is empty when no pass is found.
"""
function select_pass_times(
    u0::Vector{Float64},
    gs::GroundStation,
    opts::OrbitPropagatorOptions;
    t_end::Float64 = 24 * 3600.0,
    dt::Float64 = 10.0,
    min_elevation_deg::Float64 = 10.0,
    per_pass::Int = 8,
    max_passes::Int = 3,
)
    grid = collect(0.0:dt:t_end)
    sol = propagate_orbit(u0, (0.0, t_end), grid, opts)
    min_el = deg2rad(min_elevation_deg)

    visible = falses(length(grid))
    for i in eachindex(grid)
        r_sat = SVector{3, Float64}(sol.u[i][1:3])
        v_sat = SVector{3, Float64}(sol.u[i][4:6])
        visible[i] = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, grid[i]).elevation > min_el
    end

    times = Float64[]
    n_passes = 0
    i = 1
    while i <= length(grid) && n_passes < max_passes
        if !visible[i]
            i += 1
            continue
        end
        j = i
        while j < length(grid) && visible[j + 1]
            j += 1
        end
        idx = unique(round.(Int, range(i, j, length = min(per_pass, j - i + 1))))
        append!(times, grid[idx])
        n_passes += 1
        i = j + 1
    end
    return times, n_passes
end