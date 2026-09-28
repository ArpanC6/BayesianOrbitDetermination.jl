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
                                       omega_earth = OMEGA_EARTH)

Compute topocentric Azimuth, Elevation, Range, and Range Rate (Doppler)
from a satellite state and a ground station.

# Arguments
- `r_sat`  : Satellite ECI position (km)
- `v_sat`  : Satellite ECI velocity (km/s)
- `gs`     : `GroundStation`
- `t_sec`  : Time from epoch (seconds)

# Returns
Named tuple `(azimuth, elevation, range, range_rate)`.

This function is generic over the number type `T` (Float64, Dual, etc.),
so it can be used inside automatic-differentiation pipelines.
"""
function eci_to_station_azel_range_doppler(
    r_sat::SVector{3, T},
    v_sat::SVector{3, T},
    gs::GroundStation,
    t_sec::Float64;
    omega_earth::Float64 = OMEGA_EARTH,
) where {T}
    r_gs_ecef = ground_station_ecef(gs)
    r_gs_eci = ecef_to_eci(r_gs_ecef, t_sec; omega_earth = omega_earth)
    v_gs_eci = cross(SVector{3, Float64}(0.0, 0.0, omega_earth), r_gs_eci)

    rho_vec = r_sat - r_gs_eci
    rho_dot_vec = v_sat - v_gs_eci

    range_km = norm(rho_vec)
    range_rate = dot(rho_vec, rho_dot_vec) / range_km

    phi = deg2rad(gs.lat_deg)
    lambda = deg2rad(gs.lon_deg) + omega_earth * t_sec

    sin_lat, cos_lat = sin(phi), cos(phi)
    sin_lon, cos_lon = sin(lambda), cos(lambda)

    R_sez = SMatrix{3, 3, Float64}([
        (sin_lat * cos_lon) (sin_lat * sin_lon) (-cos_lat)
        (-sin_lon)          (cos_lon)           (0.0)
        (cos_lat * cos_lon) (cos_lat * sin_lon) (sin_lat)
    ])

    rho_sez = R_sez * rho_vec
    S, E, Z = rho_sez[1], rho_sez[2], rho_sez[3]

    elevation = asin(clamp(Z / range_km, -1.0, 1.0))
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
- `u0`               : Initial satellite state `[x, y, z, vx, vy, vz]` (km, km/s)
- `t_obs`            : Observation times (seconds from epoch)
- `gs`               : Ground station
- `sigma_range`      : Range measurement noise standard deviation (km)
- `sigma_doppler`    : Doppler measurement noise standard deviation (km/s)
- `opts`             : Propagation options
- `min_elevation_deg`: Minimum elevation for a valid pass (degrees)

# Returns
Named tuple `(t_obs, range, doppler, sol)`.

Observations are filtered to those where the satellite is above the
horizon (`elevation > min_elevation_deg`). If no observations pass the
filter, the full series is returned so that the model still has data
to fit.
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
    tspan = (t_obs[1], t_obs[end])
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

    if isempty(valid_t)
        # Fallback: use full series so the model still has observations.
        obs_r = [
            eci_to_station_azel_range_doppler(
                SVector{3, Float64}(sol.u[i][1:3]),
                SVector{3, Float64}(sol.u[i][4:6]),
                gs, t_obs[i],
            ).range + randn() * sigma_range
            for i in 1:length(t_obs)
        ]
        obs_d = [
            eci_to_station_azel_range_doppler(
                SVector{3, Float64}(sol.u[i][1:3]),
                SVector{3, Float64}(sol.u[i][4:6]),
                gs, t_obs[i],
            ).range_rate + randn() * sigma_doppler
            for i in 1:length(t_obs)
        ]
        return (t_obs = t_obs, range = obs_r, doppler = obs_d, sol = sol)
    end

    return (t_obs = valid_t, range = obs_range,
            doppler = obs_doppler, sol = sol)
end