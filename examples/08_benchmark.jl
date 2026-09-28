using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf
using Random

function run_benchmark()
    Random.seed!(42)
    d2b = Float64[]; d2m = Float64[]
    mb = Float64[]; mm = Float64[]
    cb = zeros(Int,6); cm = zeros(Int,6)
    nc = 0

    gs = GroundStation("B", 45.0, 57.5, 0.1)
    t_obs = collect(range(0.0, stop=5400.0, length=50))
    opts = OrbitPropagatorOptions(include_j2=true)
    ps = [1.0, 1.0, 1.0, 0.002, 0.002, 0.002]

    for i in 1:100
        if i % 10 == 1
            println("Trial ", i, "/100")
        end
        try
            kt = KeplerianState(6878.137+randn()*5.0, 0.001, deg2rad(51.6), deg2rad(30.0), deg2rad(40.0), deg2rad(0.0))
            ct = kepler_to_cartesian(kt)
            ut = [ct.r...; ct.v...]
            d = simulate_tracking_data(ut, t_obs, gs, 0.010, 1e-4; opts=opts)
            if length(d.t_obs) < 10
                continue
            end
            g = ut + randn(6) .* [0.5,0.5,0.5,0.001,0.001,0.001]

            bl = fit_batch_least_squares(d.t_obs, d.range, d.doppler, gs, g, 0.010, 1e-4, opts)
            push!(d2m, compute_mahalanobis_distance(ut, bl.u_estimated, bl.covariance))
            push!(mm, mean(abs.(ut[1:3] .- bl.u_estimated[1:3])))
            sv = sqrt.(max.(diag(bl.covariance), 0.0))
            for j in 1:6
                if ut[j] >= bl.u_estimated[j]-1.645*sv[j] && ut[j] <= bl.u_estimated[j]+1.645*sv[j]
                    cm[j] += 1
                end
            end

            ch = fit_bayesian_od(d.t_obs, d.range, d.doppler, gs, g, ps, opts; n_samples=500, n_adapt=250)
            zs = hcat([vec(ch[k]) for k in [:z1,:z2,:z3,:z4,:z5,:z6]]...)
            us = zs .* ps' .+ g'
            ev = EvaluateCovarianceRealism(ut, us, cov(us))
            push!(d2b, ev.mahalanobis_d2)
            push!(mb, ev.mae_position)
            cb .+= ev.coverage_90
            nc += 1
        catch e
            println("  ERROR at trial ", i, ": ", typeof(e))
        end
    end

    println("\n===== 100-TRIAL BENCHMARK =====")
    println("Completed: ", nc)
    if nc > 0
        @printf("Freq D2: %.2f +/- %.2f\n", mean(d2m), std(d2m)/sqrt(nc))
        @printf("Bayes D2: %.2f +/- %.2f\n", mean(d2b), std(d2b)/sqrt(nc))
        @printf("Freq MAE: %.4f km\n", mean(mm))
        @printf("Bayes MAE: %.4f km\n", mean(mb))
        println("Freq Coverage: ", round.((cm./nc).*100, digits=1))
        println("Bayes Coverage: ", round.((cb./nc).*100, digits=1))
        @printf("BLS overconfidence: %.0fx\n", mean(d2m)/6.0)
    end
end

run_benchmark()