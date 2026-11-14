# Frequentist / Batch Least Squares Optimizer (MLE Baseline)

"""
    fit_batch_least_squares(t_obs, range_obs, doppler_obs, gs, initial_guess,
                            sigma_range, sigma_doppler, opts;
                            max_iterations = 50, tol = 1e-9,
                            use_doppler = true)

Frequentist nonlinear batch least squares, solved with Levenberg-Marquardt.
Covariance is the Gauss-Newton approximation `cov = (J'J)^-1` of the
noise-normalized residuals, which is positive-definite by construction —
unlike the raw inverse Hessian, which is indefinite for sparse passes.

`initial_guess` is the state at epoch `t = 0`; observation times are seconds
from that epoch. Set `use_doppler = false` for range-only data such as
satellite laser ranging; `doppler_obs` is then ignored and its noise scale is
never used.

Returns a named tuple `(u_estimated, covariance, converged, iterations,
rms_residual, newton_decrement)`. `rms_residual` is the RMS of the
noise-normalized post-fit residuals and should be close to 1 when the stated
noise is correct. `converged` is true when the optimizer stopped on its own
criterion rather than exhausting `max_iterations`. `newton_decrement` is
reported for diagnostics only; it is not used to decide convergence because
the smallest eigenvalue of `J'J` is at the rounding level for sparse arcs.
"""
function fit_batch_least_squares(t_obs::Vector{Float64}, range_obs::Vector{Float64},
        doppler_obs::Vector{Float64}, gs::GroundStation,
        initial_guess::Vector{Float64}, sigma_range::Float64,
        sigma_doppler::Float64, opts::OrbitPropagatorOptions;
        max_iterations::Int = 50, tol::Float64 = 1e-9,
        use_doppler::Bool = true)

    length(t_obs) == length(range_obs) == length(doppler_obs) ||
        throw(DimensionMismatch("t_obs, range_obs and doppler_obs must have equal length"))

    function residuals(u0)
        tspan = (0.0, t_obs[end])
        sol = propagate_orbit(u0, tspan, t_obs, opts)
        T = eltype(u0)
        n = length(t_obs)
        res = Vector{T}(undef, use_doppler ? 2n : n)
        for i in eachindex(t_obs)
            r_sat = SVector{3}(sol.u[i][1], sol.u[i][2], sol.u[i][3])
            v_sat = SVector{3}(sol.u[i][4], sol.u[i][5], sol.u[i][6])
            meas = eci_to_station_azel_range_doppler(r_sat, v_sat, gs, t_obs[i])
            if use_doppler
                res[2i - 1] = (range_obs[i] - meas.range) / sigma_range
                res[2i]     = (doppler_obs[i] - meas.range_rate) / sigma_doppler
            else
                res[i] = (range_obs[i] - meas.range) / sigma_range
            end
        end
        return res
    end

    u = copy(initial_guess)
    r = residuals(u)
    cost = 0.5 * sum(abs2, r)
    damping = 1e-3
    iterations = 0
    step_norm = Inf

    for k in 1:max_iterations
        iterations = k
        J = ForwardDiff.jacobian(residuals, u)
        g = J' * r
        H = J' * J

        improved = false
        for _ in 1:12
            A = H + damping * Diagonal(diag(H)) + 1e-12 * I
            delta = -(A \ g)
            u_new = u + delta
            r_new = residuals(u_new)
            cost_new = 0.5 * sum(abs2, r_new)
            if isfinite(cost_new) && cost_new < cost
                u, r, cost = u_new, r_new, cost_new
                step_norm = norm(delta)
                damping = max(damping / 10, 1e-12)
                improved = true
                break
            end
            damping *= 10
        end

        (!improved || step_norm < tol) && break
    end

    J = ForwardDiff.jacobian(residuals, u)
    H = J' * J
    # Inverting H + eps*I directly is unstable when H is nearly singular, which
    # is normal for short arcs; the eigendecomposition keeps the covariance PSD.
    F = eigen(Symmetric(H))
    cov_est = Symmetric(F.vectors * Diagonal(1.0 ./ max.(F.values, 1e-12)) * F.vectors')
    # Newton decrement g' H^-1 g is reported as a diagnostic, but it is not used
    # to decide convergence: for sparse arcs the smallest eigenvalue of H can be
    # ~1e-10, and dividing a rounding-level gradient by it produces a meaningless
    # number. Convergence is the optimizer's own stopping rule: no damped
    # Levenberg-Marquardt step reduced the cost, or the step fell below `tol`.
    g_final = J' * r
    H_reg = F.vectors * Diagonal(max.(F.values, 1e-12)) * F.vectors'
    decrement = dot(g_final, H_reg \ g_final)
    converged = iterations < max_iterations
    rms_residual = sqrt(mean(abs2, r))

    return (u_estimated = u, covariance = cov_est, converged = converged,
            iterations = iterations, rms_residual = rms_residual,
            newton_decrement = decrement)
end
