# CCSDS Standard Conjunction Data Message (CDM) Exporter

"""
    export_ccsds_cdm(conj::ConjunctionEvent, pc::Float64, filename::String)

Export space object encounter assessment data into standardized CCSDS Conjunction Data Message (CDM) format.
"""
function export_ccsds_cdm(conj::ConjunctionEvent, pc::Float64, filename::String)
    open(filename, "w") do io
        println(io, "CCSDS_CDM_VERS = 1.0")
        println(io, "CREATION_DATE  = 2026-09-27T00:00:00.000")
        println(io, "ORIGINATOR     = BayesianOrbitDetermination.jl")
        println(io, "MESSAGE_ID     = CDM-$(conj.sat1_name)-VS-$(conj.sat2_name)")
        println(io, "")
        println(io, "TCA            = $(conj.t_tca_sec) s")
        println(io, "MISS_DISTANCE  = $(norm(conj.r2_eci - conj.r1_eci)) km")
        println(io, "COLLISION_PROB = $(pc)")
        println(io, "HARD_BODY_R    = $(conj.hard_body_radius_m) m")
        println(io, "")
        println(io, "OBJECT_1       = $(conj.sat1_name)")
        println(io, "POS_1_ECI      = [$(conj.r1_eci[1]), $(conj.r1_eci[2]), $(conj.r1_eci[3])]")
        println(io, "VEL_1_ECI      = [$(conj.v1_eci[1]), $(conj.v1_eci[2]), $(conj.v1_eci[3])]")
        println(io, "")
        println(io, "OBJECT_2       = $(conj.sat2_name)")
        println(io, "POS_2_ECI      = [$(conj.r2_eci[1]), $(conj.r2_eci[2]), $(conj.r2_eci[3])]")
        println(io, "VEL_2_ECI      = [$(conj.v2_eci[1]), $(conj.v2_eci[2]), $(conj.v2_eci[3])]")
        println(io, "=======================================================")
    end
    return filename
end
