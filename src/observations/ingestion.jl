# Observational Ingestion Module for CCSDS Tracking & TLE Ephemeris

struct GroundStationObservation
    epoch_sec::Float64
    range_km::Float64
    doppler_kms::Float64
    azimuth_rad::Float64
    elevation_rad::Float64
end

"""
    parse_ccsds_oem(oem_data::String) -> Vector{GroundStationObservation}

Parse CCSDS Orbit Ephemeris Message (OEM) formatted tracking records into structured observation array.
"""
function parse_ccsds_oem(oem_lines::Vector{String})
    observations = GroundStationObservation[]
    for line in oem_lines
        if startswith(line, "COMMENT") || startswith(line, "META") || isempty(strip(line))
            continue
        end
        tokens = split(strip(line))
        if length(tokens) >= 5
            t = parse(Float64, tokens[1])
            r = parse(Float64, tokens[2])
            d = parse(Float64, tokens[3])
            az = parse(Float64, tokens[4])
            el = parse(Float64, tokens[5])
            push!(observations, GroundStationObservation(t, r, d, az, el))
        end
    end
    return observations
end

struct CPFPrediction
    mjd::Vector{Float64}
    r_itrf_km::Matrix{Float64}
end

"""
    parse_cpf(path) -> CPFPrediction

Parse an ILRS Consolidated Prediction Format (CPF) file. Record type 10 lines
contain `MJD`, seconds of day, and ITRF position in metres; the parser returns
absolute MJD and position in km.
"""
function parse_cpf(path::String)
    mjd = Float64[]
    xyz = Float64[]
    for line in eachline(path)
        tokens = split(line)
        length(tokens) >= 8 || continue
        tokens[1] == "10" || continue
        m = parse(Float64, tokens[3])
        s = parse(Float64, tokens[4])
        push!(mjd, m + s / 86400.0)
        append!(xyz, [parse(Float64, tokens[6]) / 1000.0,
                      parse(Float64, tokens[7]) / 1000.0,
                      parse(Float64, tokens[8]) / 1000.0])
    end
    isempty(mjd) && error("no CPF record type 10 lines in $path")
    issorted(mjd) || error("CPF epochs are not monotonic in $path")
    return CPFPrediction(mjd, permutedims(reshape(xyz, 3, length(mjd))))
end

"""
    parse_slr_normal_points(paths; station = nothing) -> (mjd, range_km)

Parse ILRS CRD v2 normal-point files. For each `11` record the two-way time of
flight is converted to a vacuum range with `range = c * tof / 2`. The session
date comes from the preceding `H4` header.

The `11` epoch is defined by its epoch-event field: 0 = ground receive,
1 = spacecraft bounce, 2 = ground transmit. The returned epochs are converted
to reception time, so they can be propagated to directly and paired with
[`eci_to_station_azel_range_doppler`](@ref), which models the down-leg light
time. One-way events (3-6) are skipped.

`station` filters on the `H2` station name (case-insensitive); `nothing` keeps
all stations.
"""
function parse_slr_normal_points(paths::Vector{String};
                                 station::Union{Nothing, String} = nothing)
    mjd_all = Float64[]
    range_all = Float64[]
    for path in paths
        current_station = ""
        current_mjd0 = NaN
        for line in eachline(path)
            tokens = split(line)
            isempty(tokens) && continue
            rec = uppercase(tokens[1])
            if rec == "H2" && length(tokens) >= 2
                current_station = uppercase(tokens[2])
            elseif rec == "H4" && length(tokens) >= 5
                y = parse(Int, tokens[3])
                m = parse(Int, tokens[4])
                d = parse(Int, tokens[5])
                current_mjd0 = date_to_mjd(y, m, d)
            elseif rec == "11" && length(tokens) >= 5
                current_station == "" && continue
                (station === nothing || current_station == uppercase(station)) || continue
                isfinite(current_mjd0) || continue
                seconds = tryparse(Float64, tokens[2])
                tof = tryparse(Float64, tokens[3])
                epoch_event = tryparse(Int, tokens[5])
                (seconds === nothing || tof === nothing || epoch_event === nothing) && continue
                epoch_event in 0:2 || continue
                to_receive = epoch_event == 0 ? 0.0 : epoch_event == 1 ? tof / 2 : tof
                push!(mjd_all, current_mjd0 + (seconds + to_receive) / 86400.0)
                push!(range_all, C_LIGHT_KM_S * tof / 2.0)
            end
        end
    end
    isempty(mjd_all) && error("no normal points found" *
        (station === nothing ? "" : " for station $station"))
    perm = sortperm(mjd_all)
    return mjd_all[perm], range_all[perm]
end

"""
    date_to_mjd(year, month, day) -> Float64

Modified Julian Date at 00:00 UTC on the given calendar date.
"""
function date_to_mjd(year::Int, month::Int, day::Int)
    y, m = year, month
    if m <= 2
        y -= 1
        m += 12
    end
    A = div(y, 100)
    B = 2 - A + div(A, 4)
    jd = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + day + B - 1524.5
    return jd - 2400000.5
end
