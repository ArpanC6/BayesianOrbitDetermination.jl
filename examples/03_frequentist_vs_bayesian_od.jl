# Tutorial 03: Frequentist (batch least squares) vs Bayesian orbit determination
#
# Both estimators see the same simulated range and Doppler data. The comparison
# that matters is calibration: does the covariance each method reports actually
# cover the truth at the stated confidence?

using BayesianOrbitDetermination
using Turing
using LinearAlgebra
using Statistics
using Printf
using Random

Random.seed!(42)

println("--- BayesianOrbitDetermination.jl Tutorial 03 ---")
println("Comparison: batch least squares vs Bayesian NUTS on identical tracking data\n")

kepler_true = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
# Two passes over one station, separated by about 1.6 hours: the sparse-tracking
# regime this package is built for. The windows are fixed so the tutorial is
# deterministic; in an operational setting they come from a visibility search.
windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)

sigma_r = 0.010 # 10 m range noise
sigma_d = 1e-4  # 1 mm/s Doppler noise

opts = OrbitPropagatorOptions(include_j2 = true, include_drag = true,
                              Cd = 2.2, area_to_mass = 0.01)
data = simulate_tracking_data(u0_true, t_obs, gs, sigma_r, sigma_d; opts = opts)
@assert length(data.t_obs) == 16 "scenario produced $(length(data.t_obs)) observations"

guess_u0 = u0_true + [0.3, -0.2, 0.1, 3e-4, -2e-4, 1e-4]

println("1. Batch least squares (Levenberg-Marquardt, Gauss-Newton covariance)...")
bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                              guess_u0, sigma_r, sigma_d, opts)
@printf("   converged = %s after %d iterations, normalized residual RMS = %.3f\n",
        bls.converged, bls.iterations, bls.rms_residual)
bls_d2 = compute_mahalanobis_distance(u0_true, bls.u_estimated, bls.covariance)

println("\n2. Bayesian NUTS...")
# The posterior covariance of a sparse tracking arc is extremely ill-conditioned
# (condition number ~1e10 here). Sampling directly in state coordinates stalls;
# whitening by the BLS covariance makes the geometry isotropic.
cov_bls = Symmetric(Matrix(bls.covariance))
prior_scale = 10.0
chain = fit_bayesian_od(data.t_obs, data.range, data.doppler, gs,
                        bls.u_estimated, sqrt.(diag(cov_bls)), opts;
                        n_samples = 500, n_adapt = 300, target_accept = 0.90,
                        whitening_covariance = Matrix(cov_bls),
                        prior_scale = prior_scale,
                        sigma_range = sigma_r, sigma_doppler = sigma_d)

L = Matrix(cholesky(cov_bls).L)
z = hcat([vec(chain[k]) for k in [:z1, :z2, :z3, :z4, :z5, :z6]]...)
u_samples = (bls.u_estimated .+ prior_scale .* (L * z'))'
u_mean = vec(mean(u_samples, dims = 1))
cov_bayes = cov(u_samples)

ev = EvaluateCovarianceRealism(u0_true, u_samples, cov_bayes)
n_div = round(Int, sum(skipmissing(vec(chain[Turing.FlexiChains.Extra(:numerical_error)]))))

println("\n---------------- COMPARISON ----------------")
@printf("BLS position error      : %.4f km | D^2 = %.2f  (ideal 6)\n",
        norm(u0_true[1:3] - bls.u_estimated[1:3]), bls_d2)
@printf("Bayes position error    : %.4f km | D^2 = %.2f  (ideal 6)\n",
        norm(u0_true[1:3] - u_mean[1:3]), ev.mahalanobis_d2)
@printf("Bayes 90%% coverage      : %s\n", string(round.(ev.coverage_90 .* 100, digits = 1)))
@printf("Bayes divergences       : %d\n", n_div)
println("--------------------------------------------")
