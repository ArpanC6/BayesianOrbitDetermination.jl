# NASA NTRS Standard Covariance Realism & Goodness-of-Fit Validation Suite

"""
    compute_mahalanobis_distance(u_true::Vector{Float64}, u_est::Vector{Float64}, Cov::Matrix{Float64}) -> Float64

Compute squared Mahalanobis distance D^2 = (u_true - u_est)^T * Cov^-1 * (u_true - u_est).
Under realistic Gaussian covariance assumptions, D^2 follows a Chi-Square distribution with n degrees of freedom.
"""
function compute_mahalanobis_distance(u_true::Vector{Float64}, u_est::Vector{Float64}, Cov::Matrix{Float64})
    diff = u_true - u_est
    inv_cov = try
        inv(Cov)
    catch
        pinv(Cov)
    end
    return dot(diff, inv_cov * diff)
end

"""
    struct EvaluateCovarianceRealism

Container holding covariance realism assessment metrics for Orbit Determination:
- `mahalanobis_d2`: Squared Mahalanobis distance
- `p_value_gof`: Chi-square goodness-of-fit p-value
- `is_overconfident`: Flag indicating if estimated covariance underestimates uncertainty
- `is_underconfident`: Flag indicating if estimated covariance inflates uncertainty
- `coverage_90`: Vector of booleans indicating if ground truth fell within 90% Bayesian credible interval
"""
struct EvaluateCovarianceRealism
    mahalanobis_d2::Float64
    p_value_gof::Float64
    is_overconfident::Bool
    is_underconfident::Bool
    coverage_90::Vector{Bool}
    mae_position::Float64
    mae_velocity::Float64
end

"""
    EvaluateCovarianceRealism(u_true, posterior_samples, estimated_cov)

Assess covariance realism using Mahalanobis distance Chi-Square test and 90% credible interval coverage.
"""
function EvaluateCovarianceRealism(u_true::Vector{Float64}, posterior_samples::Matrix{Float64}, estimated_cov::Matrix{Float64})
    u_mean = vec(mean(posterior_samples, dims=1))
    d2 = compute_mahalanobis_distance(u_true, u_mean, estimated_cov)
    
    dim = length(u_true)
    chi2_dist = Chisq(dim)
    p_val = 1.0 - cdf(chi2_dist, d2)

    # 90% credible interval bounds per component
    q05 = [quantile(posterior_samples[:, j], 0.05) for j in 1:dim]
    q95 = [quantile(posterior_samples[:, j], 0.95) for j in 1:dim]
    
    coverage = [(u_true[j] >= q05[j] && u_true[j] <= q95[j]) for j in 1:dim]

    # Overconfidence / underconfidence test based on chi2 threshold (0.01 / 0.99 quantiles)
    is_over = d2 > quantile(chi2_dist, 0.99)
    is_under = d2 < quantile(chi2_dist, 0.01)

    mae_pos = mean(abs.(u_true[1:3] .- u_mean[1:3]))
    mae_vel = mean(abs.(u_true[4:6] .- u_mean[4:6]))

    return EvaluateCovarianceRealism(d2, p_val, is_over, is_under, coverage, mae_pos, mae_vel)
end
