# Extended Kalman Filter (EKF) and Unscented Kalman Filter (UKF) for Real-Time Orbit Determination

struct FilterResult
    state_history::Matrix{Float64}
    covariance_history::Vector{Matrix{Float64}}
end

"""
    run_ekf_orbit_determination(t_obs, range_obs, doppler_obs, gs, u0_prior, P0, Q, R, opts)

Run Extended Kalman Filter (EKF) sequential real-time state estimation with numerical Jacobian linearization.
"""
function run_ekf_orbit_determination(t_obs::Vector{Float64}, range_obs::Vector{Float64}, doppler_obs::Vector{Float64}, gs::GroundStation, u0_prior::Vector{Float64}, P0::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, opts::OrbitPropagatorOptions)
    n_steps = length(t_obs)
    x_est = copy(u0_prior)
    P_est = copy(P0)

    state_hist = Matrix{Float64}(undef, n_steps, 6)
    cov_hist   = Vector{Matrix{Float64}}(undef, n_steps)

    for k in 1:n_steps
        t_curr = t_obs[k]
        if k > 1
            t_prev = t_obs[k-1]
            # Predict state via numerical propagation
            sol = propagate_orbit(x_est, (t_prev, t_curr), [t_prev, t_curr], opts)
            x_pred = sol.u[end]

            # Numerical State Transition Matrix (STM) F = I + df/dx * dt
            dt = t_curr - t_prev
            f_grad(u) = orbit_dynamics(u, (opts, MU_EARTH), t_curr)
            A = ForwardDiff.jacobian(f_grad, x_pred)
            F = I(6) + A * dt

            P_pred = F * P_est * F' + Q
            x_est, P_est = x_pred, P_pred
        end

        # Measurement model prediction h(x)
        h_func(u) = begin
            r_sat = SVector{3, Float64}(u[1:3])
            v_sat = SVector{3, Float64}(u[4:6])
            m = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t_curr)
            [m.range, m.range_rate]
        end

        z_pred = h_func(x_est)
        z_meas = [range_obs[k], doppler_obs[k]]
        y = z_meas - z_pred # Innovation

        H = ForwardDiff.jacobian(h_func, x_est)
        S = H * P_est * H' + R
        K = P_est * H' * inv(S) # Kalman Gain

        x_est = x_est + K * y
        P_est = (I(6) - K * H) * P_est

        state_hist[k, :] = x_est
        cov_hist[k] = P_est
    end

    return FilterResult(state_hist, cov_hist)
end

"""
    run_ukf_orbit_determination(t_obs, range_obs, doppler_obs, gs, u0_prior, P0, Q, R, opts)

Run Unscented Kalman Filter (UKF) for non-linear orbit uncertainty propagation without Jacobian approximation.
"""
function run_ukf_orbit_determination(t_obs::Vector{Float64}, range_obs::Vector{Float64}, doppler_obs::Vector{Float64}, gs::GroundStation, u0_prior::Vector{Float64}, P0::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, opts::OrbitPropagatorOptions)
    n_dim = 6
    alpha, beta, kappa = 1e-3, 2.0, 0.0
    lambda = alpha^2 * (n_dim + kappa) - n_dim
    gamma = sqrt(n_dim + lambda)

    w_m = zeros(2n_dim + 1)
    w_c = zeros(2n_dim + 1)
    w_m[1] = lambda / (n_dim + lambda)
    w_c[1] = w_m[1] + (1 - alpha^2 + beta)
    for i in 2:(2n_dim + 1)
        w_m[i] = 0.5 / (n_dim + lambda)
        w_c[i] = w_m[i]
    end

    n_steps = length(t_obs)
    x_est = copy(u0_prior)
    P_est = copy(P0)

    state_hist = Matrix{Float64}(undef, n_steps, 6)
    cov_hist   = Vector{Matrix{Float64}}(undef, n_steps)

    for k in 1:n_steps
        t_curr = t_obs[k]
        if k > 1
            t_prev = t_obs[k-1]
            # Sigma point generation
            sqrt_P = try
                cholesky(Hermitian(P_est)).L
            catch
                sqrt(abs.(P_est))
            end

            sigma_pts = Matrix{Float64}(undef, 6, 2n_dim + 1)
            sigma_pts[:, 1] = x_est
            for i in 1:n_dim
                sigma_pts[:, i+1]       = x_est + gamma * sqrt_P[:, i]
                sigma_pts[:, i+n_dim+1] = x_est - gamma * sqrt_P[:, i]
            end

            # Propagate Sigma points
            sigma_pts_pred = Matrix{Float64}(undef, 6, 2n_dim + 1)
            for i in 1:(2n_dim + 1)
                sol = propagate_orbit(sigma_pts[:, i], (t_prev, t_curr), [t_prev, t_curr], opts)
                sigma_pts_pred[:, i] = sol.u[end]
            end

            # Predicted mean and covariance
            x_pred = zeros(6)
            for i in 1:(2n_dim + 1)
                x_pred += w_m[i] * sigma_pts_pred[:, i]
            end

            P_pred = copy(Q)
            for i in 1:(2n_dim + 1)
                diff = sigma_pts_pred[:, i] - x_pred
                P_pred += w_c[i] * (diff * diff')
            end

            x_est, P_est = x_pred, P_pred
        end

        # Measurement update
        meas = [range_obs[k], doppler_obs[k]]
        r_sat = SVector{3, Float64}(x_est[1:3])
        v_sat = SVector{3, Float64}(x_est[4:6])
        pred_meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t_curr)
        z_pred = [pred_meas.range, pred_meas.range_rate]

        H = ForwardDiff.jacobian(u -> begin
            r = SVector{3, Float64}(u[1:3]); v = SVector{3, Float64}(u[4:6])
            m = eci_to_station_azel_range_doppler(r, v, gs, t_curr)
            [m.range, m.range_rate]
        end, x_est)

        S = H * P_est * H' + R
        K = P_est * H' * inv(S)

        x_est += K * (meas - z_pred)
        P_est = (I(6) - K * H) * P_est

        state_hist[k, :] = x_est
        cov_hist[k] = P_est
    end

    return FilterResult(state_hist, cov_hist)
end
