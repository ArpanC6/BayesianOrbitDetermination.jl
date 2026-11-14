# CCSDS 502.0-B-2 Orbit Ephemeris Message (OEM) Standard Reader & Exporter

"""
    CCSDSOEMEphemeris

Structured representation of a CCSDS Orbit Ephemeris Message (OEM).
"""
struct CCSDSOEMEphemeris
    object_name::String
    object_id::String
    center_name::String
    ref_frame::String
    time_system::String
    epochs_utc::Vector{String}
    epochs_sec::Vector{Float64}
    states::Matrix{Float64} # N x 6 [x, y, z, vx, vy, vz] in km and km/s
end

"""
    export_ccsds_oem(path::String, ephem::CCSDSOEMEphemeris)

Write orbital ephemeris to disk conforming to the CCSDS 502.0-B-2 OEM standard format.
"""
function export_ccsds_oem(path::String, ephem::CCSDSOEMEphemeris)
    n_pts = size(ephem.states, 1)
    t_start = isempty(ephem.epochs_utc) ? "2026-01-01T00:00:00.000" : ephem.epochs_utc[1]
    t_stop  = isempty(ephem.epochs_utc) ? "2026-01-01T00:00:00.000" : ephem.epochs_utc[end]

    open(path, "w") do io
        println(io, "CCSDS_OEM_VERS = 2.0")
        println(io, "CREATION_DATE  = 2026-10-29T12:00:00.000")
        println(io, "ORIGINATOR     = BayesianOrbitDetermination.jl")
        println(io, "")
        println(io, "META_START")
        println(io, "OBJECT_NAME          = $(ephem.object_name)")
        println(io, "OBJECT_ID            = $(ephem.object_id)")
        println(io, "CENTER_NAME          = $(ephem.center_name)")
        println(io, "REF_FRAME            = $(ephem.ref_frame)")
        println(io, "TIME_SYSTEM          = $(ephem.time_system)")
        println(io, "START_TIME           = $(t_start)")
        println(io, "STOP_TIME            = $(t_stop)")
        println(io, "META_STOP")
        println(io, "")
        for i in 1:n_pts
            ep_str = i <= length(ephem.epochs_utc) ? ephem.epochs_utc[i] : @sprintf("%.3f", ephem.epochs_sec[i])
            @printf(io, "%-24s %16.6f %16.6f %16.6f %12.6f %12.6f %12.6f\n",
                    ep_str,
                    ephem.states[i, 1], ephem.states[i, 2], ephem.states[i, 3],
                    ephem.states[i, 4], ephem.states[i, 5], ephem.states[i, 6])
        end
    end
    return path
end

"""
    parse_ccsds_oem_file(path::String) -> CCSDSOEMEphemeris

Parse a valid CCSDS 502.0-B-2 Orbit Ephemeris Message (OEM) file containing Cartesian states.
"""
function parse_ccsds_oem_file(path::String)
    object_name = "UNKNOWN"
    object_id = "UNKNOWN"
    center_name = "EARTH"
    ref_frame = "EME2000"
    time_system = "UTC"

    epochs_utc = String[]
    epochs_sec = Float64[]
    states = Float64[]

    in_meta = false
    in_data = false

    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "COMMENT") && continue

        if occursin("META_START", s)
            in_meta = true
            continue
        elseif occursin("META_STOP", s)
            in_meta = false
            in_data = true
            continue
        end

        if in_meta
            if occursin("OBJECT_NAME", s)
                object_name = strip(split(s, "=")[end])
            elseif occursin("OBJECT_ID", s)
                object_id = strip(split(s, "=")[end])
            elseif occursin("CENTER_NAME", s)
                center_name = strip(split(s, "=")[end])
            elseif occursin("REF_FRAME", s)
                ref_frame = strip(split(s, "=")[end])
            elseif occursin("TIME_SYSTEM", s)
                time_system = strip(split(s, "=")[end])
            end
        elseif in_data
            tokens = split(s)
            if length(tokens) >= 7
                push!(epochs_utc, tokens[1])
                t_val = tryparse(Float64, tokens[1])
                push!(epochs_sec, t_val === nothing ? length(epochs_sec) * 60.0 : t_val)
                x = parse(Float64, tokens[2])
                y = parse(Float64, tokens[3])
                z = parse(Float64, tokens[4])
                vx = parse(Float64, tokens[5])
                vy = parse(Float64, tokens[6])
                vz = parse(Float64, tokens[7])
                append!(states, [x, y, z, vx, vy, vz])
            end
        end
    end

    n_records = length(epochs_utc)
    states_mat = n_records > 0 ? permutedims(reshape(states, 6, n_records)) : zeros(0, 6)
    return CCSDSOEMEphemeris(object_name, object_id, center_name, ref_frame, time_system,
                             epochs_utc, epochs_sec, states_mat)
end
