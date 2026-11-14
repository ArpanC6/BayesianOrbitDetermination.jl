# Tutorial 01: LEO orbit determination with J2, from sparse two-pass tracking.

using BayesianOrbitDetermination
using LinearAlgebra
using Statistics
using Printf
using Random

Random.seed!(42)

println()
println("BayesianOrbitDetermination.jl Tutorial 01")
println("Target: LEO state estimation from two ground-station passes")
println()

kepler_true = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
cart_true = kepler_to_cartesian(kepler_true)
u0_true = [cart_true.r...; cart_true.v...]

# Two visibility windows over one station, about 1.6 hours apart. Fixed so the
# tutorial is deterministic; operationally these come from a pass search.
gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)
opts = OrbitPropagatorOptions(include_j2 = true)

data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts = opts)
println("Generated $(length(data.t_obs)) range and Doppler observations.")

guess_u0 = u0_true + [0.3, -0.2, 0.1, 3e-4, -2e-4, 1e-4]
bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                              guess_u0, 0.010, 1e-4, opts)
@printf("BLS: converged = %s, iterations = %d, residual RMS = %.3f, position error = %.4f km\n",
        bls.converged, bls.iterations, bls.rms_residual,
        norm(u0_true[1:3] - bls.u_estimated[1:3]))

# The posterior for a sparse arc is ill-conditioned; sample in coordinates
# whitened by the BLS covariance instead of in the state itself.
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

println()
println("--- ESTIMATION RESULTS ---")
@printf("Position MAE (km)        : %.5f\n", ev.mae_position)
@printf("Velocity MAE (km/s)      : %.7f\n", ev.mae_velocity)
@printf("Mahalanobis D^2          : %.3f (chi-square p = %.4f, ideal D^2 = 6)\n",
        ev.mahalanobis_d2, ev.p_value_gof)
println("90% coverage [x, y, z, vx, vy, vz]: ", ev.coverage_90)
println("Tutorial 01 complete")
