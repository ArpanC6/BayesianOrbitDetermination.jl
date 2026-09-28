# Tutorial 01: Low Earth Orbit (LEO) Orbit Determination with J2 Perturbation

using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf

println()
println("BayesianOrbitDetermination.jl Tutorial 01")
println("Target: LEO Satellite State Estimation from Ground Station Tracking Pass")
println()

# 1. Define True LEO Satellite Initial State (Cartesian ECI)
kepler_true = KeplerianState(6878.137, 0.001, deg2rad(51.6), deg2rad(30.0), deg2rad(40.0), deg2rad(0.0))
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

@printf("True Position (km) : [%.3f, %.3f, %.3f]\n", u0_true[1], u0_true[2], u0_true[3])
@printf("True Velocity (km/s): [%.5f, %.5f, %.5f]\n\n", u0_true[4], u0_true[5], u0_true[6])

# 2. Setup Ground Station underneath satellite ground track (Overhead pass)
# Latitude 45.0° < inclination 51.6° so the satellite will pass overhead.
gs = GroundStation("LEOTrackingStation", 45.0, 57.5, 0.1)

# 3. Simulate Tracking Observations (Range & Doppler over ~100 min)
t_obs = collect(range(0.0, stop=6000.0, length=100))
sigma_range_true = 0.010   # 10 m range noise
sigma_doppler_true = 1e-4  # 1 mm/s Doppler noise

opts = OrbitPropagatorOptions(include_j2=true, include_drag=false)
data = simulate_tracking_data(u0_true, t_obs, gs, sigma_range_true, sigma_doppler_true; opts=opts)

println("Generated $(length(data.t_obs)) line-of-sight range and Doppler tracking observations.")

# 4. Setup Orbit Tracking Priors (50 m position uncertainty, 0.5 mm/s velocity uncertainty)
prior_mean = u0_true + [0.05, -0.05, 0.03, 0.0001, -0.0001, 0.0001]
prior_std  = [0.05, 0.05, 0.05, 0.0005, 0.0005, 0.0005]

# 5. Perform Bayesian Orbit Determination using Turing MCMC
println()
println("Running NUTS MCMC Bayesian Orbit Determination (target_accept=0.90)...")
chain = fit_bayesian_od(
    data.t_obs, data.range, data.doppler, gs, prior_mean, prior_std, opts;
    n_samples=1000, n_adapt=500, target_accept=0.90,
)

# 6. Extract Posterior Estimates & Covariance Realism Validation
z_samples = hcat([vec(chain[k]) for k in [:z1, :z2, :z3, :z4, :z5, :z6]]...)
u_samples = z_samples .* prior_std' .+ prior_mean'

u_est_mean = vec(mean(u_samples, dims=1))
cov_est = cov(u_samples)

eval_realism = EvaluateCovarianceRealism(u0_true, u_samples, cov_est)

println()
println("--- ESTIMATION RESULTS ---")
@printf("Position MAE (km)  : %.5f\n", eval_realism.mae_position)
@printf("Velocity MAE (km/s) : %.7f\n", eval_realism.mae_velocity)
@printf("Mahalanobis Distance D^2: %.3f (Chi-Square GOF p-val: %.4f)\n",
        eval_realism.mahalanobis_d2, eval_realism.p_value_gof)
println("90% Credible Interval Coverage for [x, y, z, vx, vy, vz]: ",
        eval_realism.coverage_90)
println()
println("Tutorial 01 complete successfully")