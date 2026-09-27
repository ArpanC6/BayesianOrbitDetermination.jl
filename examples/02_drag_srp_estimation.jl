# Tutorial 02: Orbit Determination with Atmospheric Drag & Solar Radiation Pressure (SRP)

using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf

println("=== BayesianOrbitDetermination.jl Tutorial 02 ===")
println("Target: Orbit Estimation under Combined Perturbations (J2 + Drag + SRP)\n")

kepler_true = KeplerianState(6778.137, 0.002, deg2rad(45.0), deg2rad(0.0), deg2rad(0.0), deg2rad(0.0))
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

gs = GroundStation("Goldstone", 35.426, -116.890, 1.0)
t_obs = collect(range(0.0, stop=5400.0, length=25))

opts = OrbitPropagatorOptions(include_j2=true, include_drag=true, include_srp=true, Cd=2.2, Cr=1.3, area_to_mass=0.015)
data = simulate_tracking_data(u0_true, t_obs, gs, 0.015, 2e-4; opts=opts)

prior_mean = u0_true + [0.5, -0.5, 0.5, 0.001, -0.001, 0.001]
prior_std  = [2.0, 2.0, 2.0, 0.003, 0.003, 0.003]

println("Sampling posterior for high-altitude LEO with Drag & SRP...")
chain = fit_bayesian_od(t_obs, data.range, data.doppler, gs, prior_mean, prior_std, opts; n_samples=300, n_adapt=100)

u_samples = hcat([vec(chain[k]) for k in [:u1, :u2, :u3, :u4, :u5, :u6]]...)
u_est_mean = vec(mean(u_samples, dims=1))
cov_est = cov(u_samples)

eval_realism = EvaluateCovarianceRealism(u0_true, u_samples, cov_est)

@printf("Position MAE: %.4f km | Velocity MAE: %.6f km/s\n", eval_realism.mae_position, eval_realism.mae_velocity)
@printf("Mahalanobis D^2: %.2f | Coverage: %s\n", eval_realism.mahalanobis_d2, string(eval_realism.coverage_90))
println("Tutorial 02 complete!")
