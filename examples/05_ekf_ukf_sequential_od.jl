# Tutorial 05: Sequential EKF and UKF Orbit Determination

using BayesianOrbitDetermination
using LinearAlgebra
using Printf
using Random

Random.seed!(42)

println("--- BayesianOrbitDetermination.jl Tutorial 05 ---")
println("Sequential EKF vs UKF over two ground-station passes\n")

# Same two-pass scenario as Tutorials 01-03.
kepler_true = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)
opts = OrbitPropagatorOptions(include_j2 = true)

data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts = opts)
println("Generated $(length(data.t_obs)) range and Doppler observations over two passes\n")

u0_prior = u0_true + [0.5, -0.5, 0.5, 5e-4, -5e-4, 5e-4]
P0 = Matrix(Diagonal([0.25, 0.25, 0.25, 1e-6, 1e-6, 1e-6]))
Q  = Matrix(Diagonal(fill(1e-14, 6)))
R  = Matrix(Diagonal([0.010^2, (1e-4)^2]))

println("Running sequential EKF (variational STM)...")
ekf_res = run_ekf_orbit_determination(data.t_obs, data.range, data.doppler, gs, u0_prior, P0, Q, R, opts)

println("Running sequential UKF (sigma points)...")
ukf_res = run_ukf_orbit_determination(data.t_obs, data.range, data.doppler, gs, u0_prior, P0, Q, R, opts)

t_final = data.t_obs[end]
u_true_final = propagate_orbit(u0_true, (0.0, t_final), [t_final], opts).u[end]

ekf_pos_err = norm(u_true_final[1:3] - ekf_res.state_history[end, 1:3])
ukf_pos_err = norm(u_true_final[1:3] - ukf_res.state_history[end, 1:3])
ekf_sigma = sqrt(sum(diag(ekf_res.covariance_history[end][1:3, 1:3])))
ukf_sigma = sqrt(sum(diag(ukf_res.covariance_history[end][1:3, 1:3])))

@printf("EKF final position error: %.4f km (3D 1-sigma %.4f km)\n", ekf_pos_err, ekf_sigma)
@printf("UKF final position error: %.4f km (3D 1-sigma %.4f km)\n", ukf_pos_err, ukf_sigma)
@printf("EKF-UKF state difference: %.2e km\n",
        norm(ekf_res.state_history[end, 1:3] - ukf_res.state_history[end, 1:3]))
println("Tutorial 05 complete")
