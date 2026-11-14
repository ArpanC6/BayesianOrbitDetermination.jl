# ISRO & Indian Space Situational Awareness (SSA) Ground Tracking Network Definitions
#
# References:
# - ISTRAC (ISRO Telemetry, Tracking and Command Network), Bengaluru
# - Project NETRA (Network for Space Object Tracking and Analysis)
# - PRL Mount Abu Satellite Laser Ranging (SLR) Observatory, Gurushikhar

"""
    const ISRO_GROUND_STATIONS

Canonical geodetic coordinates (lat deg, lon deg, alt km) for ISRO tracking, telemetry,
laser ranging, and SSA radar stations across India and overseas downrange facilities.
"""
const ISRO_GROUND_STATIONS = Dict{Symbol, GroundStation}(
    :BANGALORE_ISTRAC   => GroundStation("ISTRAC_Bengaluru", 12.9716, 77.5946, 0.920),
    :SRIHARIKOTA_SHAR   => GroundStation("SDSC_SHAR_Sriharikota", 13.7199, 80.2307, 0.020),
    :PORT_BLAIR_ISTRAC  => GroundStation("ISTRAC_Port_Blair", 11.6234, 92.7265, 0.030),
    :LUCKNOW_ISTRAC     => GroundStation("ISTRAC_Lucknow", 26.8467, 80.9462, 0.120),
    :MOUNT_ABU_SLR      => GroundStation("PRL_Mount_Abu_SLR", 24.5926, 72.7156, 1.680),
    :PONMUDI_NETRA      => GroundStation("NETRA_Radar_Ponmudi", 8.7600, 77.1100, 0.950),
    :MAURITIUS_ISTRAC   => GroundStation("ISTRAC_Mauritius", -20.2444, 57.7556, 0.050),
    :BRUNEI_ISTRAC      => GroundStation("ISTRAC_Brunei", 4.9031, 114.9398, 0.025),
    :BHARTI_ANTARCTICA  => GroundStation("ISRO_Bharti_Antarctica", -69.4069, 76.1872, 0.035)
)

"""
    get_isro_station(station_id::Symbol) -> GroundStation

Retrieve predefined ISRO ground station by symbolic identifier, e.g. `:BANGALORE_ISTRAC` or `:SRIHARIKOTA_SHAR`.
"""
function get_isro_station(station_id::Symbol)
    haskey(ISRO_GROUND_STATIONS, station_id) || error("Unknown ISRO station :$station_id. Available: $(keys(ISRO_GROUND_STATIONS))")
    return ISRO_GROUND_STATIONS[station_id]
end
