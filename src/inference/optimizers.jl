# Frequentist / Batch Least Squares Optimizer (MLE Baseline)

"""
    fit_batch_least_squares(t_obs, range_obs, doppler_obs, gs, initial_guess, sigma_range, sigma_doppler, opts)

Frequentist Nonlinear Batch Least Squares (MAP / Maximum Likelihood Estimation) for baseline comparison against standard orbit determination tools (e.g., Orekit / MONTE).
"""
function fit_batch_least_squares(t_obs::Vector{Float64}, range_obs::Vector{Float64}, doppler_obs::Vector{Float64}, gs::GroundStation, initial_guess::Vector{Float64}, sigma_range::Float64, sigma_doppler::Float64, opts::OrbitPropagatorOptions)
    
    function loss(u0)
        tspan = (t_obs[1], t_obs[end])
        sol = propagate_orbit(u0, tspan, t_obs, opts)

        n = length(t_obs)
        res = 0.0
        for i in 1:n
            r_sat = SVector{3}(sol.u[i][1:3]) # Fixed SVector
            v_sat = SVector{3}(sol.u[i][4:6]) # Fixed SVector
            meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t_obs[i])

            r_res = (range_obs[i] - meas.range) / sigma_range
            d_res = (doppler_obs[i] - meas.range_rate) / sigma_doppler
            res += r_res^2 + d_res^2
        end
        return 0.5 * res
    end

    res = optimize(loss, initial_guess, LBFGS(), Optim.Options(iterations=200, g_tol=1e-6))
    u_est = Optim.minimizer(res)

    # Approximate covariance from Hessian at MLE
    H = ForwardDiff.hessian(loss, u_est)
    cov_est = try
        inv(H)
    catch
        pinv(H)
    end

    return (u_estimated=u_est, covariance=cov_est, result=res)
end