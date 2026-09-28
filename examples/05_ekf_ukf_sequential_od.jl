# Tutorial 05: Real-Time Extended Kalman Filter (EKF) vs Unscented Kalman Filter (UKF)

using BayesianOrbitDetermination
using LinearAlgebra
using Printf

println("=== BayesianOrbitDetermination.jl Tutorial 05 ===")
println("Comparison: Sequential EKF vs UKF Real-Time State Estimation\n")

kepler_true = KeplerianState(6878.137, 0.001, deg2rad(51.6), deg2rad(30.0), deg2rad(40.0), deg2rad(0.0))
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

gs = GroundStation("LEOTrackingStation", 45.0, 57.5, 0.1)
t_obs = collect(range(0.0, stop=5400.0, length=100))

opts = OrbitPropagatorOptions(include_j2=true)
data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts=opts)

n_obs = length(data.t_obs)
println("Generated $n_obs observations")

u0_prior = u0_true + [0.1, -0.1, 0.05, 0.0002, -0.0002, 0.0001]
P0 = Matrix(Diagonal([0.5, 0.5, 0.5, 1e-4, 1e-4, 1e-4]))
Q  = Matrix(Diagonal([1e-4, 1e-4, 1e-4, 1e-6, 1e-6, 1e-6]))
R  = Matrix(Diagonal([0.010^2, (1e-4)^2]))

println("Running Sequential EKF...")
ekf_res = run_ekf_orbit_determination(data.t_obs, data.range, data.doppler, gs, u0_prior, P0, Q, R, opts)

println("Running Sequential UKF...")
ukf_res = run_ukf_orbit_determination(data.t_obs, data.range, data.doppler, gs, u0_prior, P0, Q, R, opts)

# Propagate true state to final observation time for correct comparison
t_final = data.t_obs[end]
sol_true = propagate_orbit(u0_true, (data.t_obs[1], t_final), [t_final], opts)
u_true_final = sol_true.u[end]

ekf_pos_err = norm(u_true_final[1:3] - ekf_res.state_history[end, 1:3])
ukf_pos_err = norm(u_true_final[1:3] - ukf_res.state_history[end, 1:3])

@printf("Final EKF Position Error: %.4f km\n", ekf_pos_err)
@printf("Final UKF Position Error: %.4f km\n", ukf_pos_err)
println("Tutorial 05 complete!")
