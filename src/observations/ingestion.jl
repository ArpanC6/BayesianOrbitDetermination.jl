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
