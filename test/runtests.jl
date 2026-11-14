using Test
using BayesianOrbitDetermination
using DelimitedFiles
using ForwardDiff
using LinearAlgebra
using StaticArrays
using Statistics
using Random

Random.seed!(42)

@testset "BayesianOrbitDetermination.jl" begin

    @testset "Keplerian <-> Cartesian Conversion" begin
        k_orig = KeplerianState(7000.0, 0.01, deg2rad(28.5), deg2rad(45.0), deg2rad(30.0), deg2rad(15.0))
        c_state = kepler_to_cartesian(k_orig)
        k_reconv = cartesian_to_kepler(c_state)

        @test isapprox(k_orig.a, k_reconv.a, rtol=1e-5)
        @test isapprox(k_orig.e, k_reconv.e, rtol=1e-5)
    end

    @testset "Physics: inverse-square gravity and closed orbit" begin
        r = SVector(7000.0, 0.0, 0.0)
        a_vec = BayesianOrbitDetermination.two_body_acceleration(r)
        @test norm(a_vec) ≈ MU_EARTH / 7000.0^2 rtol = 1e-14

        v_circ = sqrt(MU_EARTH / 7000.0)
        u0 = [7000.0, 0.0, 0.0, 0.0, v_circ, 0.0]
        period = 2pi * sqrt(7000.0^3 / MU_EARTH)

        sol = propagate_orbit(u0, (0.0, period), [period],
                              OrbitPropagatorOptions(include_j2 = false))
        @test sol.u[end] ≈ u0 rtol = 1e-8
    end

    @testset "Physics: energy conservation" begin
        u0 = [7000.0, 0.0, 0.0, 0.0, 7.0, 1.0]
        energy(u) = 0.5 * sum(abs2, u[4:6]) - MU_EARTH / norm(u[1:3])
        t = collect(range(0.0, 3600.0, length = 25))

        sol = propagate_orbit(u0, (0.0, 3600.0), t,
                              OrbitPropagatorOptions(include_j2 = false))
        E0 = energy(u0)
        @test all(abs(energy(u) - E0) / abs(E0) < 1e-8 for u in sol.u)
    end

    @testset "Numerical propagation matches the analytic Kepler solution" begin
        kep = KeplerianState(7000.0, 0.001, deg2rad(51.6), deg2rad(30.0),
                             deg2rad(40.0), deg2rad(10.0))
        cart = kepler_to_cartesian(kep)
        u0 = [cart.r...; cart.v...]

        # Mean anomaly at epoch, then invert Kepler's equation at each time.
        n = sqrt(MU_EARTH / kep.a^3)
        E0 = 2 * atan(sqrt((1 - kep.e) / (1 + kep.e)) * tan(kep.nu / 2))
        M0 = E0 - kep.e * sin(E0)

        opts = OrbitPropagatorOptions(include_j2 = false, reltol = 1e-12, abstol = 1e-12)
        for dt in (3600.0, 6 * 3600.0, 24 * 3600.0)
            M = M0 + n * dt
            E = M
            for _ in 1:30
                E -= (E - kep.e * sin(E) - M) / (1 - kep.e * cos(E))
            end
            nu = 2 * atan(sqrt((1 + kep.e) / (1 - kep.e)) * tan(E / 2))
            analytic = kepler_to_cartesian(
                KeplerianState(kep.a, kep.e, kep.i, kep.RAAN, kep.argp, nu))

            sol = propagate_orbit(u0, (0.0, dt), [dt], opts)
            err_km = norm(SVector{3}(sol.u[end][1:3]...) - analytic.r)
            @test err_km < 1e-3   # 1 m over up to one day; measured ~2e-6 km
        end
    end

    @testset "J2 secular nodal drift matches the analytic rate" begin
        kep = KeplerianState(7000.0, 0.001, deg2rad(51.6), deg2rad(30.0),
                             deg2rad(40.0), 0.0)
        cart = kepler_to_cartesian(kep)
        u0 = [cart.r...; cart.v...]
        opts = OrbitPropagatorOptions(include_j2 = true, reltol = 1e-11, abstol = 1e-11)

        mean_motion = sqrt(MU_EARTH / kep.a^3)
        p = kep.a * (1 - kep.e^2)
        rate = -1.5 * mean_motion * J2_EARTH * (R_EARTH / p)^2 * cos(kep.i)
        t_end = 50 * 2pi / mean_motion

        sol = propagate_orbit(u0, (0.0, t_end), [0.0, t_end], opts)
        k0 = cartesian_to_kepler(CartesianState(
            SVector{3}(sol.u[1][1:3]...), SVector{3}(sol.u[1][4:6]...)))
        k1 = cartesian_to_kepler(CartesianState(
            SVector{3}(sol.u[2][1:3]...), SVector{3}(sol.u[2][4:6]...)))
        dnode = mod(k1.RAAN - k0.RAAN + pi, 2pi) - pi

        @test isapprox(dnode, rate * t_end, rtol = 0.01)
    end

    @testset "Zonal harmonics: J2 limit and potential gradient" begin
        r = SVector(3205.83, 4997.61, 3488.02)
        @test zonal_harmonic_perturbation(r; n_max = 2) ≈ j2_perturbation(r) rtol = 1e-12

        # a = -grad U with U = mu/r * sum J_n (R/r)^n P_n(z/r)
        function potential(x)
            rn = norm(x)
            u = x[3] / rn
            s = 0.0
            P_prev2, P_prev1 = 1.0, u
            all_j = (J2_EARTH, J3_EARTH, J4_EARTH, J5_EARTH)
            for n in 2:5
                Pn = n == 2 ? (3u^2 - 1) / 2 :
                     ((2n - 1) * u * P_prev1 - (n - 1) * P_prev2) / n
                s += all_j[n - 1] * (R_EARTH / rn)^n * Pn
                P_prev2, P_prev1 = P_prev1, Pn
            end
            return MU_EARTH / rn * s
        end
        @test zonal_harmonic_perturbation(r; n_max = 5) ≈
              -ForwardDiff.gradient(potential, collect(r)) rtol = 1e-10
    end

    @testset "TLE Parsing & SGP4 Dynamics" begin
        l1 = "1 25544U 98067A   26270.50000000  .00016717  00000-0  30000-3 0  9993"
        l2 = "2 25544  51.6400 208.1234 0004500  65.1234 295.0000 15.49500000420002"
        tle = parse_tle(l1, l2, "ISS")
        @test tle.norad_id == 25544
        @test isapprox(tle.inc_deg, 51.64, rtol=1e-3)

        c_prop = propagate_sgp4(tle, 600.0)
        @test norm(c_prop.r) > 6000.0

        # Regression vector from SatelliteToolboxSgp4 3.0.0 (Vallado SGP4);
        # this is the real SGP4 output in TEME, not a Kepler + J2 approximation.
        c30 = propagate_sgp4(tle, 1800.0)
        @test isapprox(c30.r, SVector(4445.342816822122, -1907.3060317047552,
                                      4762.608537289984); atol = 1e-6)
        @test isapprox(c30.v, SVector(5.0647079290803925, 5.0935045710169184,
                                      -2.6769198389541216); atol = 1e-9)

        k0 = tle_to_keplerian(tle)
        @test 6700.0 < k0.a < 6900.0
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

        ekf_res = run_ekf_orbit_determination(data.t_obs, data.range, data.doppler, gs, u0_true, P0, Q, R, opts)
        @test size(ekf_res.state_history, 1) == length(data.t_obs)
        @test size(ekf_res.covariance_history, 1) == length(data.t_obs)
        @test eigmin(Symmetric(ekf_res.covariance_history[end])) > -1e-12

        ukf_res = run_ukf_orbit_determination(data.t_obs, data.range, data.doppler, gs, u0_true, P0, Q, R, opts)
        @test size(ukf_res.state_history, 1) == length(data.t_obs)
    end

    @testset "Variational state-transition matrix" begin
        u0 = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
        opts_fd = OrbitPropagatorOptions(include_j2 = true, reltol = 1e-12, abstol = 1e-12)
        _, Phi = propagate_state_and_stm(u0, (0.0, 600.0), opts_fd)

        for j in 1:6
            delta = j <= 3 ? 1e-3 : 1e-4
            u_plus = copy(u0); u_plus[j] += delta
            p_plus = propagate_orbit(u_plus, (0.0, 600.0), [600.0], opts_fd).u[end]
            p_base = propagate_orbit(u0, (0.0, 600.0), [600.0], opts_fd).u[end]
            @test isapprox(Phi[:, j], (p_plus .- p_base) ./ delta, rtol = 1e-3, atol = 1e-5)
        end
    end

    @testset "EKF covariance prediction on a two-pass arc" begin
        kepler = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
        cart = kepler_to_cartesian(kepler)
        u0 = [cart.r...; cart.v...]
        gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
        opts = OrbitPropagatorOptions(include_j2 = true)
        windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
        t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)
        data = simulate_tracking_data(u0, t_obs, gs, 0.010, 1e-4; opts = opts)

        P0 = Matrix(Diagonal([0.25, 0.25, 0.25, 1e-6, 1e-6, 1e-6]))
        Q  = Matrix(Diagonal(fill(1e-14, 6)))
        R  = Matrix(Diagonal([0.01^2, (1e-4)^2]))
        x0 = u0 .+ [0.5, -0.5, 0.5, 5e-4, -5e-4, 5e-4]

        res = run_ekf_orbit_determination(data.t_obs, data.range, data.doppler, gs, x0, P0, Q, R, opts)
        truth_end = propagate_orbit(u0, (0.0, data.t_obs[end]), [data.t_obs[end]], opts).u[end]
        pos_err = norm(res.state_history[end, 1:3] - truth_end[1:3])
        pos_sigma = sqrt(sum(diag(res.covariance_history[end][1:3, 1:3])))

        @test pos_err < 0.1
        @test pos_err < 3 * pos_sigma
        @test pos_sigma < 0.2
    end

    @testset "Conjunction Assessment Pc & CDM Export" begin
        r1 = SVector(6800.0, 0.0, 0.0); v1 = SVector(0.0, 7.6, 0.0)
        r2 = SVector(6800.1, 0.02, 0.01); v2 = SVector(0.0, -7.6, 0.5)
        cov1 = Matrix(Diagonal(fill(0.01, 6)))
        cov2 = Matrix(Diagonal(fill(0.01, 6)))
        conj = ConjunctionEvent("SAT1", "SAT2", 500.0, r1, v1, cov1, r2, v2, cov2, 10.0)

        pc = compute_collision_probability_foster(conj)
        @test pc >= 0.0 && pc <= 1.0

        tmp_summary = joinpath(tempdir(), "test_conjunction_summary.txt")
        export_cdm_summary(conj, pc, tmp_summary)
        @test isfile(tmp_summary)
        @test occursin("NOT a CCSDS", read(tmp_summary, String))

        # Full CCSDS 508.0-B-1 CDM export
        tmp_cdm = joinpath(tempdir(), "test_conjunction.cdm")
        export_ccsds_cdm(conj, pc, tmp_cdm)
        @test isfile(tmp_cdm)
        cdm_content = read(tmp_cdm, String)
        @test occursin("CCSDS_CDM_VERS = 1.0", cdm_content)
        @test occursin("OBJECT = OBJECT1", cdm_content)
        @test occursin("OBJECT = OBJECT2", cdm_content)
        @test occursin("COLLISION_PROBABILITY", cdm_content)

        # Monte Carlo Collision Probability validation
        pc_mc, pc_mc_se = compute_collision_probability_monte_carlo(conj; n_samples = 5000)
        @test 0.0 <= pc_mc <= 1.0
        @test isapprox(pc, pc_mc, atol = 3 * max(pc_mc_se, 0.02))
    end

    @testset "CCSDS 502.0-B-2 OEM Ephemeris Export & Parse" begin
        epochs = ["2026-10-29T12:00:00.000", "2026-10-29T12:01:00.000"]
        t_sec = [0.0, 60.0]
        states = [7000.0 0.0 0.0 0.0 7.5 0.0;
                  7000.5 10.0 0.0 0.0 7.4 0.0]
        oem = CCSDSOEMEphemeris("TESTSAT", "2026-999A", "EARTH", "EME2000", "UTC", epochs, t_sec, states)
        tmp_oem = joinpath(tempdir(), "test_orbit.oem")
        export_ccsds_oem(tmp_oem, oem)
        @test isfile(tmp_oem)

        parsed = parse_ccsds_oem_file(tmp_oem)
        @test parsed.object_name == "TESTSAT"
        @test parsed.object_id == "2026-999A"
        @test size(parsed.states) == (2, 6)
        @test parsed.states[1, 1] ≈ 7000.0
    end

    @testset "ISRO Ground Stations & Earth Orientation" begin
        shar = get_isro_station(:SRIHARIKOTA_SHAR)
        @test shar.name == "SDSC_SHAR_Sriharikota"
        @test isapprox(shar.lat_deg, 13.72, atol = 0.1)

        bengaluru = get_isro_station(:BANGALORE_ISTRAC)
        @test bengaluru.name == "ISTRAC_Bengaluru"

        # GMST test
        gmst = gmst_rad(0.0)
        @test 0.0 <= gmst <= 2pi

        # ITRF to ECI conversion
        r_itrf = SVector(6378.137, 0.0, 0.0)
        r_eci = itrf_to_eci_state(r_itrf, 0.0)
        @test isapprox(norm(r_eci), 6378.137, atol = 1e-4)
    end

    @testset "NASA LaRC UQ Challenge Solver" begin
        u_true = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        samples = randn(50, 6) .* 0.01 .+ u_true'
        res = solve_nasa_larc_uq_subproblem(samples, u_true)
        @test res.reliability_index >= 0.0
    end

    @testset "Tracking data use u0 at epoch t = 0" begin
        u0 = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
        gs = GroundStation("Equator", 0.0, 0.0, 0.0)
        opts = OrbitPropagatorOptions(include_j2 = true)
        t_obs = collect(range(30.0, 330.0, length = 10))

        data = simulate_tracking_data(u0, t_obs, gs, 0.0, 0.0; opts = opts)
        @test length(data.t_obs) >= 5

        sol = propagate_orbit(u0, (0.0, maximum(data.t_obs)), data.t_obs, opts)
        for k in eachindex(data.t_obs)
            r = SVector{3}(sol.u[k][1], sol.u[k][2], sol.u[k][3])
            v = SVector{3}(sol.u[k][4], sol.u[k][5], sol.u[k][6])
            expected = eci_to_station_azel_range_doppler(r, v, gs, data.t_obs[k])
            @test data.range[k] ≈ expected.range atol = 1e-9
            @test data.doppler[k] ≈ expected.range_rate atol = 1e-12
        end
    end

    @testset "ILRS CRD and CPF parsing" begin
        @test date_to_mjd(2000, 1, 1) == 51544.0
        @test date_to_mjd(2025, 6, 9) == 60835.0
        @test date_to_mjd(2024, 2, 29) == 60369.0

        data_dir = joinpath(@__DIR__, "..", "validation", "slr", "data")
        cpf = parse_cpf(joinpath(data_dir, "starlette_cpf_250609_16101.sgf"))
        @test cpf.mjd[1] == 60835.0
        @test issorted(cpf.mjd)
        @test size(cpf.r_itrf_km, 2) == 3
        @test 7300.0 < norm(cpf.r_itrf_km[1, :]) < 7600.0  # Starlette orbit radius

        mjd, rng = parse_slr_normal_points(
            [joinpath(data_dir, "starlette_20250609.np2")]; station = "YARL")
        @test length(mjd) > 20
        @test issorted(mjd)
        @test all(500.0 .< rng .< 4500.0)
        @test 60835.0 < mjd[1] < 60836.0

        # CRD record 11 epoch events: 0 = receive, 1 = bounce, 2 = transmit.
        # The parser returns reception epochs.
        mktempdir() do dir
            path = joinpath(dir, "epochs.np2")
            open(path, "w") do io
                println(io, "h2 TEST       7090  5 13 3       ILRS")
                println(io, "h4 1 2025 6 9 0 0 0 2025 6 9 1 0 0 180 1")
                for event in (0, 1, 2, 3)
                    println(io, "11  1000.000000000000 0.010000000000 cfg1 $(event)   30.0   10  40.0  0.1 -0.7 0.9 1.0 0 0.0")
                end
            end
            mjd_ev, rng_ev = parse_slr_normal_points([path]; station = "TEST")
            @test length(mjd_ev) == 3          # one-way event 3 skipped
            @test rng_ev ≈ fill(BayesianOrbitDetermination.C_LIGHT_KM_S * 0.005, 3)
            @test (mjd_ev[1] - 60835.0) * 86400 ≈ 1000.0           # receive
            @test (mjd_ev[2] - 60835.0) * 86400 ≈ 1000.005         # bounce + tof/2
            @test (mjd_ev[3] - 60835.0) * 86400 ≈ 1000.010         # transmit + tof
        end
    end

    @testset "Batch least squares recovers the truth" begin
        k = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
        c = kepler_to_cartesian(k)
        u0 = [c.r...; c.v...]
        gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
        opts = OrbitPropagatorOptions(include_j2 = true, include_drag = true,
                                      Cd = 2.2, area_to_mass = 0.01)
        windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
        t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)

        data = simulate_tracking_data(u0, t_obs, gs, 0.010, 1e-4; opts = opts)
        @test length(data.t_obs) == 16

        guess = u0 + [0.3, -0.2, 0.1, 3e-4, -2e-4, 1e-4]
        bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                                      guess, 0.010, 1e-4, opts)
        @test bls.converged
        @test bls.rms_residual < 1.5
        @test norm(bls.u_estimated[1:3] - u0[1:3]) < 0.5
    end

    @testset "Covariance-whitened Bayesian OD" begin
        k = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
        c = kepler_to_cartesian(k)
        u0 = [c.r...; c.v...]
        gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
        opts = OrbitPropagatorOptions(include_j2 = true, include_drag = true,
                                      Cd = 2.2, area_to_mass = 0.01)
        windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
        t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)

        data = simulate_tracking_data(u0, t_obs, gs, 0.010, 1e-4; opts = opts)
        guess = u0 + [0.3, -0.2, 0.1, 3e-4, -2e-4, 1e-4]
        bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                                      guess, 0.010, 1e-4, opts)

        C = Symmetric(Matrix(bls.covariance))
        L = Matrix(cholesky(C).L)
        chain = fit_bayesian_od(data.t_obs, data.range, data.doppler, gs,
                                bls.u_estimated, sqrt.(diag(C)), opts;
                                n_samples = 100, n_adapt = 100, target_accept = 0.90,
                                whitening_covariance = Matrix(C), prior_scale = 10.0)

        z = hcat([vec(chain[n]) for n in (:z1, :z2, :z3, :z4, :z5, :z6)]...)
        u = (bls.u_estimated .+ 10.0 .* (L * z'))'
        @test size(u) == (100, 6)
        @test all(isfinite, u)

        ev = EvaluateCovarianceRealism(u0, u, cov(u))
        @test isfinite(ev.mahalanobis_d2)
        @test ev.mahalanobis_d2 >= 0.0
        @test norm(u0[1:3] - vec(mean(u, dims = 1))[1:3]) < 0.5
    end

    @testset "Covariance realism diagnostics" begin
        Random.seed!(11)
        C = Matrix(Diagonal([1.0, 2.0, 3.0, 0.1, 0.1, 0.1]))
        L = cholesky(Symmetric(C)).L

        # A calibrated covariance puts D^2 at chi-square(6) on average.
        d2 = [compute_mahalanobis_distance(L * randn(6), zeros(6), C) for _ in 1:200]
        @test 5.2 < mean(d2) < 6.8

        truth = L * randn(6)
        samples = (L * randn(6, 4000))'
        ev = EvaluateCovarianceRealism(truth, samples, C)
        @test count(ev.coverage_90) >= 4
        @test 0.0 <= ev.p_value_gof <= 1.0

        # Ill-conditioned covariances must not produce a negative squared distance.
        C_ill = Matrix(Diagonal([1e-18, 1e-18, 1e-18, 1.0, 1.0, 1.0]))
        @test compute_mahalanobis_distance(ones(6), zeros(6), C_ill) >= 0.0
    end

    @testset "SBC uniformity test" begin
        Random.seed!(3)
        n = 1000
        n_draws = 500

        # Ranks live on 0:n_draws; binning must use the draw count, not the
        # number of simulations.
        uniform_ranks = rand(0:n_draws, n, 6)
        res = sbc_uniformity_test(uniform_ranks; n_draws = n_draws)
        @test all(res.p_value .> 0.001)
        @test res.dof == 9

        # Ranks confined to the lower half of 0:n_draws are not uniform.
        clustered = rand(0:(n_draws ÷ 2), n, 6)
        res_clustered = sbc_uniformity_test(clustered; n_draws = n_draws)
        @test all(res_clustered.p_value .< 1e-6)

        degenerate_ranks = zeros(Int, n, 6)
        res_bad = sbc_uniformity_test(degenerate_ranks)
        @test all(res_bad.p_value .< 1e-6)
    end

    @testset "SBC run end to end" begin
        u0 = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
        gs = GroundStation("Equator", 0.0, 0.0, 0.0)
        opts = OrbitPropagatorOptions(include_j2 = true)
        prior_std = [0.1, 0.1, 0.1, 1e-4, 1e-4, 1e-4]
        prior_cov = Matrix(Diagonal(prior_std .^ 2))

        res = run_sbc(gs, opts, u0, prior_cov;
                      n_simulations = 1, prior_scale = 1.0,
                      infer_noise = true, n_samples = 40, n_adapt = 40,
                      t_end = 6000.0, dt = 30.0, per_pass = 6,
                      max_passes = 2, seed = 7, verbose = false)
        @test size(res.ranks) == (1, 8)
        @test res.param_names == [:z1, :z2, :z3, :z4, :z5, :z6,
                                  :log_sigma_range, :log_sigma_doppler]
        @test all(0 .<= res.ranks .<= 40)
        @test all(isfinite, res.rhat)
        @test all(isfinite, res.ess)
    end

    @testset "Independent Orekit cross-validation" begin
        u0 = [3205.830222953533, 4997.609661816432, 3488.022461182947,
              -6.035351552351067, 0.6629103392316869, 4.597259758177123]
        ref_dir = joinpath(@__DIR__, "..", "validation", "orekit", "reference")
        for (file, drag) in (("reference_j2.csv", false),
                             ("reference_j2drag.csv", true))
            data = readdlm(joinpath(ref_dir, file), ',', Float64; skipstart = 2)
            t = data[:, 1]
            ref = data[:, 2:7]
            opts = OrbitPropagatorOptions(include_j2 = true, include_drag = drag,
                                          reltol = 1e-12, abstol = 1e-12)
            sol = propagate_orbit(u0, (0.0, t[end]), collect(t), opts)
            got = reduce(vcat, permutedims.(collect.(sol.u)))
            dp = maximum(norm(got[i, 1:3] .- ref[i, 1:3]) for i in eachindex(t))
            dv = maximum(norm(got[i, 4:6] .- ref[i, 4:6]) for i in eachindex(t))
            # Orekit DP8(5,3) reference; measured differences are a few mm.
            @test dp < 0.05  # m
            @test dv < 1e-5  # m/s
        end
    end

end
