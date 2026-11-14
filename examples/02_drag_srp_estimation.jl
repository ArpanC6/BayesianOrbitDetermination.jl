# Tutorial 02: LEO orbit determination with J2, drag, and solar radiation pressure.

using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf
using Random

Random.seed!(42)

println("--- BayesianOrbitDetermination.jl Tutorial 02 ---")
println("Target: estimation under J2 + drag + SRP from two passes\n")

kepler_true = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)

opts = OrbitPropagatorOptions(include_j2 = true, include_drag = true,
                              include_srp = true, Cd = 2.2, Cr = 1.3,
                              area_to_mass = 0.015)
data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts = opts)
println("Generated $(length(data.t_obs)) observations.")

guess_u0 = u0_true + [0.3, -0.2, 0.1, 3e-4, -2e-4, 1e-4]
bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                              guess_u0, 0.010, 1e-4, opts)
@printf("BLS: converged = %s, residual RMS = %.3f, position error = %.4f km\n",
        bls.converged, bls.rms_residual, norm(u0_true[1:3] - bls.u_estimated[1:3]))

C = Symmetric(Matrix(bls.covariance))
L = Matrix(cholesky(C).L)
prior_scale = 10.0
chain = fit_bayesian_od(data.t_obs, data.range, data.doppler, gs,
                        bls.u_estimated, sqrt.(diag(C)), opts;
                        n_samples = 500, n_adapt = 300, target_accept = 0.90,
                        whitening_covariance = Matrix(C), prior_scale = prior_scale)

z = hcat([vec(chain[k]) for k in [:z1, :z2, :z3, :z4, :z5, :z6]]...)
u_samples = (bls.u_estimated .+ prior_scale .* (L * z'))'
ev = EvaluateCovarianceRealism(u0_true, u_samples, cov(u_samples))

@printf("\nPosition MAE: %.4f km | Velocity MAE: %.6f km/s\n",
        ev.mae_position, ev.mae_velocity)
@printf("Mahalanobis D^2: %.2f (ideal 6) | Coverage: %s\n",
        ev.mahalanobis_d2, string(ev.coverage_90))
println("Tutorial 02 complete")
