using Test
using BayesianOrbitDetermination
using LinearAlgebra
using StaticArrays
using Statistics

@testset "BayesianOrbitDetermination.jl Full 12-Month Test Suite" begin

    @testset "Keplerian <-> Cartesian Conversion" begin
        k_orig = KeplerianState(7000.0, 0.01, deg2rad(28.5), deg2rad(45.0), deg2rad(30.0), deg2rad(15.0))
        c_state = kepler_to_cartesian(k_orig)
        k_reconv = cartesian_to_kepler(c_state)

        @test isapprox(k_orig.a, k_reconv.a, rtol=1e-5)
        @test isapprox(k_orig.e, k_reconv.e, rtol=1e-5)
    end

    @testset "TLE Parsing & SGP4 Dynamics" begin
        l1 = "1 25544U 98067A   26270.50000000  .00016717  00000-0  30000-3 0  9993"
        l2 = "2 25544  51.6400 208.1234 0004500  65.1234 295.0000 15.49500000420002"
        tle = parse_tle(l1, l2, "ISS")
        @test tle.norad_id == 25544
        @test isapprox(tle.inc_deg, 51.64, rtol=1e-3)

        c_prop = propagate_sgp4(tle, 600.0)
        @test norm(c_prop.r) > 6000.0
    end

    @testset "Orbit Dynamics & Propagation" begin
        u0 = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
        tspan = (0.0, 3600.0)
        t_save = [0.0, 1800.0, 3600.0]
        opts = OrbitPropagatorOptions(include_j2=true)

        sol = propagate_orbit(u0, tspan, t_save, opts)
        @test length(sol.u) == 3
    end

    @testset "Sequential EKF / UKF Filters" begin
        u0_true = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
        gs = GroundStation("TestStation", 0.0, 0.0, 0.0)
        t_obs = collect(range(0.0, stop=600.0, length=5))
        opts = OrbitPropagatorOptions(include_j2=true)
        data = simulate_tracking_data(u0_true, t_obs, gs, 0.01, 1e-4; opts=opts)

        P0 = Matrix(Diagonal(fill(1.0, 6)))
        Q  = Matrix(Diagonal(fill(1e-6, 6)))
        R  = Matrix(Diagonal([0.01^2, (1e-4)^2]))

        ekf_res = run_ekf_orbit_determination(t_obs, data.range, data.doppler, gs, u0_true, P0, Q, R, opts)
        @test size(ekf_res.state_hist, 1) == 5
    end

    @testset "Conjunction Assessment Pc & CDM Export" begin
        r1 = SVector(6800.0, 0.0, 0.0); v1 = SVector(0.0, 7.6, 0.0)
        r2 = SVector(6800.1, 0.02, 0.01); v2 = SVector(0.0, -7.6, 0.5)
        cov1 = Matrix(Diagonal(fill(0.01, 6)))
        cov2 = Matrix(Diagonal(fill(0.01, 6)))
        conj = ConjunctionEvent("SAT1", "SAT2", 500.0, r1, v1, cov1, r2, v2, cov2, 10.0)

        pc = compute_collision_probability_foster(conj)
        @test pc >= 0.0 && pc <= 1.0

        tmp_cdm = joinpath(tempdir(), "test_cdm.cdm")
        export_ccsds_cdm(conj, pc, tmp_cdm)
        @test isfile(tmp_cdm)
    end

    @testset "NASA LaRC UQ Challenge Solver" begin
        u_true = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        samples = randn(50, 6) .* 0.01 .+ u_true'
        res = solve_nasa_larc_uq_subproblem(samples, u_true)
        @test res.reliability_index >= 0.0
    end

end
