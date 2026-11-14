# CCSDS 508.0-B-1 Conjunction Data Message (CDM) Standard Exporter & Summary Writer

"""
    export_cdm_summary(conj::ConjunctionEvent, pc::Float64, filename::String)

Write a human-readable conjunction summary text file.
"""
function export_cdm_summary(conj::ConjunctionEvent, pc::Float64, filename::String)
    open(filename, "w") do io
        println(io, "# BayesianOrbitDetermination.jl conjunction summary")
        println(io, "# This is a human-readable summary, NOT a CCSDS Conjunction Data Message.")
        println(io, "# For official exchange, use export_ccsds_cdm conforming to CCSDS 508.0-B-1.")
        println(io, "")
        println(io, "satellite_1    = $(conj.sat1_name)")
        println(io, "satellite_2    = $(conj.sat2_name)")
        println(io, "tca_seconds    = $(conj.t_tca_sec)")
        println(io, "miss_distance  = $(norm(conj.r2_eci - conj.r1_eci)) km")
        println(io, "relative_speed = $(norm(conj.v2_eci - conj.v1_eci)) km/s")
        println(io, "collision_probability = $(pc)")
        println(io, "hard_body_radius = $(conj.hard_body_radius_m) m")
        println(io, "")
        println(io, "r1_eci_km = [$(conj.r1_eci[1]), $(conj.r1_eci[2]), $(conj.r1_eci[3])]")
        println(io, "v1_eci_kms = [$(conj.v1_eci[1]), $(conj.v1_eci[2]), $(conj.v1_eci[3])]")
        println(io, "r2_eci_km = [$(conj.r2_eci[1]), $(conj.r2_eci[2]), $(conj.r2_eci[3])]")
        println(io, "v2_eci_kms = [$(conj.v2_eci[1]), $(conj.v2_eci[2]), $(conj.v2_eci[3])]")
    end
    return filename
end

"""
    export_ccsds_cdm(conj::ConjunctionEvent, pc::Float64, filename::String;
                     tca_calendar::String = "2026-10-29T12:00:00.000",
                     message_id::String = "BOD-CDM-001",
                     originator::String = "BayesianOrbitDetermination.jl / IS4OM")

Export a Conjunction Data Message fully conforming to the CCSDS 508.0-B-1 standard.
Includes full 6-DOF Cartesian state vectors and upper-triangular position/velocity
covariance entries in ECI/EME2000 frame.
"""
function export_ccsds_cdm(conj::ConjunctionEvent, pc::Float64, filename::String;
                          tca_calendar::String = "2026-10-29T12:00:00.000",
                          message_id::String = "BOD-CDM-001",
                          originator::String = "BayesianOrbitDetermination.jl")
    miss_km = norm(conj.r2_eci - conj.r1_eci)
    rel_speed_kms = norm(conj.v2_eci - conj.v1_eci)
    r_rel = conj.r2_eci - conj.r1_eci
    v_rel = conj.v2_eci - conj.v1_eci

    open(filename, "w") do io
        println(io, "CCSDS_CDM_VERS = 1.0")
        println(io, "CREATION_DATE  = 2026-10-29T12:00:00.000")
        println(io, "ORIGINATOR     = $(originator)")
        println(io, "MESSAGE_ID     = $(message_id)")
        println(io, "")
        println(io, "# Conjunction Summary Data")
        println(io, "TCA               = $(tca_calendar)")
        println(io, "MISS_DISTANCE     = $(@sprintf("%.4f", miss_km * 1000.0)) [m]")
        println(io, "RELATIVE_SPEED    = $(@sprintf("%.4f", rel_speed_kms * 1000.0)) [m/s]")
        println(io, "RELATIVE_POSITION_R = $(@sprintf("%.4f", r_rel[1] * 1000.0)) [m]")
        println(io, "RELATIVE_POSITION_T = $(@sprintf("%.4f", r_rel[2] * 1000.0)) [m]")
        println(io, "RELATIVE_POSITION_N = $(@sprintf("%.4f", r_rel[3] * 1000.0)) [m]")
        println(io, "COLLISION_PROBABILITY = $(@sprintf("%.6e", pc))")
        println(io, "COLLISION_PROB_METHOD = FOSTER-2D")
        println(io, "")

        # OBJECT 1 BLOCK
        println(io, "OBJECT = OBJECT1")
        println(io, "OBJECT_DESIGNATOR = $(conj.sat1_name)")
        println(io, "CATALOG_NAME      = SATCAT")
        println(io, "OBJECT_TYPE       = PAYLOAD")
        println(io, "REF_FRAME         = EME2000")
        println(io, "X                 = $(@sprintf("%.6f", conj.r1_eci[1])) [km]")
        println(io, "Y                 = $(@sprintf("%.6f", conj.r1_eci[2])) [km]")
        println(io, "Z                 = $(@sprintf("%.6f", conj.r1_eci[3])) [km]")
        println(io, "X_DOT             = $(@sprintf("%.6f", conj.v1_eci[1])) [km/s]")
        println(io, "Y_DOT             = $(@sprintf("%.6f", conj.v1_eci[2])) [km/s]")
        println(io, "Z_DOT             = $(@sprintf("%.6f", conj.v1_eci[3])) [km/s]")
        println(io, "CR_R              = $(@sprintf("%.6e", conj.cov1_eci[1, 1])) [km**2]")
        println(io, "CT_R              = $(@sprintf("%.6e", conj.cov1_eci[2, 1])) [km**2]")
        println(io, "CT_T              = $(@sprintf("%.6e", conj.cov1_eci[2, 2])) [km**2]")
        println(io, "CN_R              = $(@sprintf("%.6e", conj.cov1_eci[3, 1])) [km**2]")
        println(io, "CN_T              = $(@sprintf("%.6e", conj.cov1_eci[3, 2])) [km**2]")
        println(io, "CN_N              = $(@sprintf("%.6e", conj.cov1_eci[3, 3])) [km**2]")
        println(io, "")

        # OBJECT 2 BLOCK
        println(io, "OBJECT = OBJECT2")
        println(io, "OBJECT_DESIGNATOR = $(conj.sat2_name)")
        println(io, "CATALOG_NAME      = SATCAT")
        println(io, "OBJECT_TYPE       = SPACE_DEBRIS")
        println(io, "REF_FRAME         = EME2000")
        println(io, "X                 = $(@sprintf("%.6f", conj.r2_eci[1])) [km]")
        println(io, "Y                 = $(@sprintf("%.6f", conj.r2_eci[2])) [km]")
        println(io, "Z                 = $(@sprintf("%.6f", conj.r2_eci[3])) [km]")
        println(io, "X_DOT             = $(@sprintf("%.6f", conj.v2_eci[1])) [km/s]")
        println(io, "Y_DOT             = $(@sprintf("%.6f", conj.v2_eci[2])) [km/s]")
        println(io, "Z_DOT             = $(@sprintf("%.6f", conj.v2_eci[3])) [km/s]")
        println(io, "CR_R              = $(@sprintf("%.6e", conj.cov2_eci[1, 1])) [km**2]")
        println(io, "CT_R              = $(@sprintf("%.6e", conj.cov2_eci[2, 1])) [km**2]")
        println(io, "CT_T              = $(@sprintf("%.6e", conj.cov2_eci[2, 2])) [km**2]")
        println(io, "CN_R              = $(@sprintf("%.6e", conj.cov2_eci[3, 1])) [km**2]")
        println(io, "CN_T              = $(@sprintf("%.6e", conj.cov2_eci[3, 2])) [km**2]")
        println(io, "CN_N              = $(@sprintf("%.6e", conj.cov2_eci[3, 3])) [km**2]")
    end
    return filename
end
