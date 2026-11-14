# Real-data validation: ILRS satellite laser ranging to Starlette (NORAD 07646).
#
# Data: EUROLAS Data Center (DGFI-TUM), CRD v2 normal points from station
# Yarragadee (YARL, 7090) on 2025-06-09..11, and the matching CPF prediction.
# See README.md for provenance and licensing.
#
# The script performs three independent checks:
#   1. CPF prediction vs the SLR normal points (external orbit check).
#   2. Our force model propagated from the CPF state vs the SLR normal points.
#   3. A range-only batch least-squares fit to the SLR normal points.
# Every number in RESULTS.md is written by this script.

using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf
using StaticArrays

const DATA = joinpath(@__DIR__, "data")
const CPF_FILE = joinpath(DATA, "starlette_cpf_250609_16101.sgf")
const NP_FILES = [joinpath(DATA, "starlette_20250609.np2"),
                  joinpath(DATA, "starlette_20250610.np2"),
                  joinpath(DATA, "starlette_20250611.np2")]
const STATION = "YARL"
const STATION_LAT = -29.0464   # deg, ILRS
const STATION_LON = 115.3467   # deg east, ILRS
const STATION_ALT = 0.244      # km
const MIN_ELEV_DEG = 15.0      # reject low-elevation data (refraction, weighting)
const STARLETTE_MASS = 47.3    # kg
const STARLETTE_RADIUS = 0.12  # m
const STARLETTE_CD = 2.1
const STARLETTE_CR = 1.1

# --- interpolation on the 180 s CPF grid -----------------------------------
function lagrange_weights(tn)
    n = length(tn)
    w = ones(n)
    for k in 1:n, j in 1:n
        j == k && continue
        w[k] /= (tn[k] - tn[j])
    end
    return w
end

function lagrange_eval(t, y, tq; m::Int = 8)
    n = length(t)
    out = zeros(length(tq))
    for (q, t0) in enumerate(tq)
        i = searchsortedfirst(t, t0)
        lo = clamp(i - m ÷ 2, 1, n - m + 1)
        idx = lo:lo + m - 1
        w = lagrange_weights(@view t[idx])
        num = 0.0
        den = 0.0
        hit = false
        for (k, ik) in enumerate(idx)
            if t0 == t[ik]
                out[q] = y[ik]
                hit = true
                break
            end
            c = w[k] / (t0 - t[ik])
            num += c * y[ik]
            den += c
        end
        hit || (out[q] = num / den)
    end
    return out
end

function lagrange_derivative_at_node(t, y, i0; m::Int = 8)
    n = length(t)
    lo = clamp(i0 - m ÷ 2, 1, n - m + 1)
    idx = collect(lo:lo + m - 1)
    tn = t[idx]
    yn = y[idx]
    w = lagrange_weights(tn)
    k0 = findfirst(==(i0), idx)
    d = zeros(m)
    for k in 1:m
        k == k0 && continue
        d[k] = (w[k] / w[k0]) / (tn[k0] - tn[k])
    end
    d[k0] = -sum(d)
    return dot(d, yn)
end

function lagrange_eval_deriv(t, y, tq; m::Int = 8)
    n = length(t)
    val = zeros(length(tq))
    der = zeros(length(tq))
    for (q, t0) in enumerate(tq)
        i = searchsortedfirst(t, t0)
        lo = clamp(i - m ÷ 2, 1, n - m + 1)
        idx = lo:lo + m - 1
        w = lagrange_weights(@view t[idx])
        hit = false
        for (k, ik) in enumerate(idx)
            if t0 == t[ik]
                val[q] = y[ik]
                der[q] = lagrange_derivative_at_node(t, y, ik; m = m)
                hit = true
                break
            end
        end
        hit && continue
        num = 0.0
        den = 0.0
        dnum = 0.0
        dden = 0.0
        for (k, ik) in enumerate(idx)
            c = w[k] / (t0 - t[ik])
            num += c * y[ik]
            den += c
            dnum -= c * y[ik] / (t0 - t[ik])
            dden -= c / (t0 - t[ik])
        end
        val[q] = num / den
        der[q] = (dnum * den - num * dden) / den^2
    end
    return val, der
end

function describe(x)
    return (n = length(x), mean_m = 1e3 * mean(x), rms_m = 1e3 * sqrt(mean(abs2, x)),
            max_m = 1e3 * maximum(abs, x))
end

function main()
    cpf = parse_cpf(CPF_FILE)
    mjd_np, range_np = parse_slr_normal_points(NP_FILES; station = STATION)
    t_cpf = (cpf.mjd .- cpf.mjd[1]) .* 86400.0
    t_np = (mjd_np .- cpf.mjd[1]) .* 86400.0

    # Order-of-accuracy check for the CPF interpolant: compare 8- and 10-point
    # evaluations at grid midpoints, where neither stencil is exact. The 8-point
    # stencil degenerates to one-sided near the ends of the file, so report the
    # interior and the edge separately.
    let mid = t_cpf[21:end-21] .+ 90.0, mid_all = t_cpf[1:end-1] .+ 90.0
        err(tq) = maximum(maximum(abs, lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], tq; m = 8) .-
                                   lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], tq; m = 10))
                          for j in 1:3)
        @printf("CPF interpolation: interior max |8-point - 10-point| = %.3f mm (%d midpoints)\n",
                1e6 * err(mid), length(mid))
        @printf("                   edge max |8-point - 10-point| = %.3f mm (one-sided stencils)\n",
                1e6 * err(mid_all))
    end

    # CPF position and velocity at every normal-point epoch, in ITRF.
    r_cpf = hcat([lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], t_np) for j in 1:3]...)
    v_cpf = hcat([lagrange_eval_deriv(t_cpf, cpf.r_itrf_km[:, j], t_np)[2] for j in 1:3]...)

    # Model range through the same two-way light-time measurement model used in
    # the fit. CRD record 11 in this data set is transmit-tagged (epoch event 2),
    # so parse_slr_normal_points has already converted the epochs to reception
    # time and the CPF is interpolated there.
    rng_cpf = zeros(length(t_np))
    elev = zeros(length(t_np))
    dh_range = zeros(length(t_np))
    for i in eachindex(t_np)
        jd_i = cpf.mjd[1] + t_np[i] / 86400.0
        r_eci, v_eci, th_i = itrf_to_eci_state(
            SVector{3}(r_cpf[i, :]), SVector{3}(v_cpf[i, :]), jd_i)
        gs_i = GroundStation("$(STATION)-inertial", STATION_LAT,
                             STATION_LON + rad2deg(th_i), STATION_ALT)
        m = eci_to_station_azel_range_doppler(r_eci, v_eci, gs_i, 0.0)
        m_geom = eci_to_station_azel_range_doppler(r_eci, v_eci, gs_i, 0.0;
                                                   light_time = false)
        rng_cpf[i] = m.range
        elev[i] = m.elevation
        dh_range[i] = m.range - m_geom.range
    end

    keep = findall(>(deg2rad(MIN_ELEV_DEG)), elev)
    t_np = t_np[keep]
    range_np = range_np[keep]
    mjd_np = mjd_np[keep]
    r_cpf = r_cpf[keep, :]
    v_cpf = v_cpf[keep, :]
    rng_cpf = rng_cpf[keep]
    elev = elev[keep]
    dh_range = dh_range[keep]

    res_cpf = range_np .- rng_cpf

    println("CPF vs SLR: ", describe(res_cpf))
    @printf("  light-time range correction: rms %.1f m, max %.1f m\n",
            1e3 * sqrt(mean(abs2, dh_range)), 1e3 * maximum(abs, dh_range))

    # Diagnostic for the epoch convention: use each transmit-tagged epoch as if
    # it were already a reception time. The residual should be much worse.
    let
        t_tx = t_np .- 2.0 .* range_np ./ BayesianOrbitDetermination.C_LIGHT_KM_S
        rr = hcat([lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], t_tx) for j in 1:3]...)
        vv = hcat([lagrange_eval_deriv(t_cpf, cpf.r_itrf_km[:, j], t_tx)[2] for j in 1:3]...)
        r_wrong = zeros(length(t_tx))
        for i in eachindex(t_tx)
            jd_i = cpf.mjd[1] + t_tx[i] / 86400.0
            r_eci, v_eci, th_i = itrf_to_eci_state(
                SVector{3}(rr[i, :]), SVector{3}(vv[i, :]), jd_i)
            gs_i = GroundStation("$(STATION)-inertial", STATION_LAT,
                                 STATION_LON + rad2deg(th_i), STATION_ALT)
            r_wrong[i] = eci_to_station_azel_range_doppler(r_eci, v_eci, gs_i, 0.0).range
        end
        @printf("  diagnostic: without the epoch-event conversion, rms %.1f m\n",
                1e3 * sqrt(mean(abs2, range_np .- r_wrong)))
    end

    # Interpolation order check at the data epochs themselves.
    let
        d = maximum(maximum(abs, lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], t_np; m = 8) .-
                              lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], t_np; m = 10))
                    for j in 1:3)
        @printf("  CPF interpolation at the %d SLR epochs: max |8-point - 10-point| = %.3f mm\n",
                length(t_np), 1e6 * d)
    end
    @printf("  elevation range: %.1f - %.1f deg\n", rad2deg(minimum(elev)), rad2deg(maximum(elev)))
    @printf("  span: %.2f days, %d normal points\n",
            (mjd_np[end] - mjd_np[1]), length(range_np))
    flush(stdout)

    # --- our force model from the CPF state --------------------------------
    i0 = searchsortedlast(t_cpf, t_np[1])
    t0 = t_cpf[i0]
    r0_itrf = SVector{3}([lagrange_eval(t_cpf, cpf.r_itrf_km[:, j], [t0])[1] for j in 1:3])
    v0_itrf = SVector{3}([lagrange_derivative_at_node(t_cpf, cpf.r_itrf_km[:, j], i0)
                          for j in 1:3])
    jd0 = cpf.mjd[1] + t0 / 86400.0
    r0_eci, v0_eci, theta0 = itrf_to_eci_state(r0_itrf, v0_itrf, jd0)
    u0 = [r0_eci...; v0_eci...]

    gs_eci = GroundStation("$(STATION)-inertial", STATION_LAT,
                           STATION_LON + rad2deg(theta0), STATION_ALT)
    area_to_mass = pi * STARLETTE_RADIUS^2 / STARLETTE_MASS
    opts = OrbitPropagatorOptions(include_j2 = true, include_j3j4 = true,
                                  include_drag = true,
                                  include_srp = true,
                                  Cd = STARLETTE_CD, Cr = STARLETTE_CR,
                                  area_to_mass = area_to_mass,
                                  reltol = 1e-11, abstol = 1e-11,
                                  epoch_jd = jd0)

    t_obs = t_np .- t0
    sol = propagate_orbit(u0, (0.0, t_obs[end]), t_obs, opts)
    range_model = [eci_to_station_azel_range_doppler(
                       SVector{3}(sol.u[i][1:3]), SVector{3}(sol.u[i][4:6]),
                       gs_eci, t_obs[i]).range for i in eachindex(t_obs)]
    res_model = range_np .- range_model
    println("Model from CPF state: ", describe(res_model))
    day = floor.(Int, mjd_np .- mjd_np[1])
    for d in unique(day)
        sel = day .== d
        r = res_model[sel]
        @printf("  day +%d: n=%2d  mean %8.1f m  rms %8.1f m  max %8.1f m\n",
                d, length(r), 1e3 * mean(r), 1e3 * sqrt(mean(abs2, r)), 1e3 * maximum(abs, r))
    end

    # Effect of each perturbation at the end of the arc, relative to two-body+J2.
    let
        base = propagate_orbit(u0, (0.0, t_obs[end]), [t_obs[end]],
                               OrbitPropagatorOptions(include_j2 = true,
                                                      reltol = 1e-11, abstol = 1e-11))
        r_base = base.u[end][1:3]
        for (label, o) in (("+J3/J4", OrbitPropagatorOptions(
                                include_j2 = true, include_j3j4 = true,
                                reltol = 1e-11, abstol = 1e-11)),
                           ("+J3/J4+drag", OrbitPropagatorOptions(
                                include_j2 = true, include_j3j4 = true,
                                include_drag = true, Cd = STARLETTE_CD,
                                area_to_mass = area_to_mass,
                                reltol = 1e-11, abstol = 1e-11)),
                           ("+J3/J4+drag+SRP", opts))
            s = propagate_orbit(u0, (0.0, t_obs[end]), [t_obs[end]], o)
            @printf("  %-16s final position difference: %.3f m\n", label,
                    1e3 * norm(s.u[end][1:3] - r_base))
        end
    end

    # --- range-only batch least squares ------------------------------------
    function fit_arc(sel, label)
        t = t_obs[sel]
        r = range_np[sel]
        bls = fit_batch_least_squares(t, r, zeros(length(t)), gs_eci,
                                      u0, 0.010, 1.0, opts;
                                      use_doppler = false, max_iterations = 25,
                                      tol = 1e-7)
        sol_fit = propagate_orbit(bls.u_estimated, (0.0, t[end]), t, opts)
        r_fit = [eci_to_station_azel_range_doppler(
                     SVector{3}(sol_fit.u[i][1:3]), SVector{3}(sol_fit.u[i][4:6]),
                     gs_eci, t[i]).range for i in eachindex(t)]
        res = r .- r_fit
        println("BLS $(label): ", describe(res))
        @printf("  converged=%s iterations=%d rms=%.3f (normalized at 10 m)\n",
                bls.converged, bls.iterations, bls.rms_residual)
        @printf("  state correction: %.3f m, %.4f mm/s\n",
                1e3 * norm(bls.u_estimated[1:3] - u0[1:3]),
                1e6 * norm(bls.u_estimated[4:6] - u0[4:6]))
        sig = sqrt.(diag(bls.covariance))
        @printf("  formal 1-sigma (10 m weight): pos %.2f m, vel %.4f mm/s\n",
                1e3 * norm(sig[1:3]), 1e6 * norm(sig[4:6]))
        # The residual RMS is the square root of the reduced chi-square in units
        # of the assumed 10 m noise; scale the covariance a posteriori by it.
        a_post = sig .* bls.rms_residual
        @printf("  a-posteriori 1-sigma (x%.2f): pos %.2f m, vel %.4f mm/s\n",
                bls.rms_residual, 1e3 * norm(a_post[1:3]), 1e6 * norm(a_post[4:6]))
        @printf("  first/last residual: %.1f / %.1f m\n", 1e3 * res[1], 1e3 * res[end])
        for d in unique(day[sel])
            s2 = day[sel] .== d
            rr = res[s2]
            @printf("    day +%d: n=%2d  mean %7.1f m  rms %7.1f m  max %7.1f m\n",
                    d, length(rr), 1e3 * mean(rr), 1e3 * sqrt(mean(abs2, rr)),
                    1e3 * maximum(abs, rr))
        end
        flush(stdout)
        return (label = label, bls = bls, res = res, n = length(t),
                rms = sqrt(mean(abs2, res)))
    end

    one_day = day .== 0
    fit_arc(one_day, "first day")
    fit_arc(trues(length(t_obs)), "full arc")
end

main()
