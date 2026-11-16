# BayesianOrbitDetermination.jl

Bayesian orbit determination, covariance realism, and simulation-based calibration for Julia and the SciML ecosystem.

[![Build Status](https://github.com/ArpanC6/BayesianOrbitDetermination.jl/workflows/CI/badge.svg)](https://github.com/ArpanC6/BayesianOrbitDetermination.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## The question this package exists to answer

Under sparse ground-station tracking, does a full Bayesian posterior give a better-calibrated state covariance than the linearized covariance from a batch least-squares fit? Most published comparisons cannot answer this, because the optimizer has not converged, or because the sampler cannot move in an ill-conditioned posterior, and the resulting covariance is then blamed on the method rather than on the computation.

This package is built so that both estimators can be computed correctly and compared against the same simulated truth: a Levenberg-Marquardt batch fit with a Gauss-Newton covariance and an explicit convergence flag, and NUTS in Turing.jl with a covariance-whitened parameterization. An earlier version of this README quoted a "150x overconfident" linearized covariance; that number came from an unconverged optimizer, and it has been withdrawn.

## Headline result

50 sparse-tracking LEO scenarios, each with up to three passes and 24 range and Doppler observations over 24 h, fitted independently by Levenberg-Marquardt batch least squares and by a whitened NUTS posterior on the same data:

| Metric | BLS | Bayesian NUTS | Ideal |
|---|---|---|---|
| Mean Mahalanobis D^2 (6-DOF) | 7.37 +/- 0.74 | 7.89 +/- 0.80 | 6 |
| Per-axis position MAE (km) | 0.0263 +/- 0.0037 | 0.0262 +/- 0.0037 | 0 |
| Per-component 90% coverage (%) | 92, 92, 94, 90, 92, 92 | 94, 90, 90, 90, 92, 90 | 90 |
| Divergent chains | -- | 0 / 50 | 0 |
| Worst r-hat / smallest ESS | -- | 1.046 / 395 | 1.0 / large |

Marginal coverage is at the nominal rate for both estimators, and the two estimators produce statistically indistinguishable position errors (26.3 m and 26.2 m per-axis MAE). The joint D^2 is 1.9 standard errors above 6 for BLS and 2.4 for the Bayesian posterior; the paired difference (Bayes minus BLS) is +0.52 +/- 0.20, so the full posterior does not remove the residual. Whether that residual is a small overconfidence of the linearized covariance in the least-constrained direction, non-Gaussianity of the posterior, or seed variance cannot be settled by 50 replications; the benchmark exists to settle it.

What the study does establish is negative. Once the optimizer is converged and the sampler is whitened, there is no detectable Bayesian advantage in the point estimate or the marginal calibration on this problem. The earlier "150x overconfident" claim was an optimizer artifact and has been withdrawn.

Reproduce it with:

    julia --project=. benchmarks/run_covariance_realism_study.jl

The script writes `benchmarks/RESULTS.md` and `benchmarks/results_per_trial.csv` with every trial's Mahalanobis distance, position error, convergence flag, residual RMS, divergence count, r-hat, and effective sample size. Nothing in the table is typed by hand.

## One-command demo

    julia --project=. run_demo.jl

This prints the validated ILRS SLR residual (7.83 m RMS), the Orekit cross-check (1.2 mm / 2.8 mm), the covariance-realism diagnostics (D^2 = 7.37 / 7.89, zero divergences), and the conjunction Pc. Every number is printed by the script, not typed by hand.

## Real tracking data: satellite laser ranging to Starlette

The benchmark above is self-consistent: the same models generate and fit the data. `validation/slr/` removes that for the measurement model by fitting real ILRS laser ranges to Starlette from Yarragadee (station 7090), 2025-06-09 to 2025-06-11, against the independent ILRS Consolidated Prediction Format ephemeris.

- The CRD normal points are transmit-tagged (epoch event 2), so the parser converts each epoch to reception time before the two-way light-time solution. Applying the range model without that conversion inflates the CPF residual from **7.83 m to 39.5 m RMS**, which is how the convention was verified.
- Through the same measurement model, the package agrees with the CPF ephemeris over 98 normal points to **7.83 m RMS** (mean 4.15 m, max 19.1 m) at 15-66 deg elevation over 2.74 days. The light-time correction itself is 20.9 m RMS.
- A range-only Levenberg-Marquardt fit with the package's reduced force model reaches 24.7 m RMS on the first day and 58.8 m RMS over the full arc.
- Starting from the CPF state, our reduced force model diverges from the CPF by about 4 km RMS over 2.74 days. That is the expected cost of omitting tesseral harmonics, solid Earth tides, relativity, and atmospheric delay, and it is why this is a measurement-model validation, not an operational POD accuracy claim.

Every number is printed by `julia --project=. validation/slr/run_starlette.jl`; `validation/slr/README.md` records the data provenance, checksums, and limitations.

## What is implemented and validated

- Exact two-body gravity, plus J2-J5 zonal harmonics (Legendre recurrence, enabled with `include_j3j4 = true`), exponential atmospheric drag, and SRP with an analytic solar ephemeris and a conical shadow. The numerical two-body propagation matches the closed-form Kepler solution to 2 mm over 24 h at `reltol = abstol = 1e-12` (about 2 m at the default `1e-9`), the J2 nodal drift matches the analytic secular rate to better than 1% over 50 orbits, the zonal acceleration is checked against a finite-difference gradient of the potential, and closed-orbit periodicity and energy conservation are tested in `test/runtests.jl`.
- The J2 and J2 + drag propagator is cross-validated against **Orekit 13.1.9** on the benchmark scenario: the maximum position difference over 24 h is 1.2 mm for J2 and 2.8 mm with the exponential atmosphere and Cd = 2.2. The Orekit-generated reference states, the Java generator, and the checksums are committed in `validation/orekit/`, and the comparison runs in CI.
- Ground-station range and range-rate (Doppler) models with WGS-84 geodetic to ECEF to ECI conversion, including the two-way light-time solution: at LEO ranges the bounce epoch is several milliseconds before reception, which is a ~20 m effect on an SLR range and is solved by fixed-point iteration. The initial state is defined at epoch `t = 0`, and a regression test checks that simulated observations match a direct propagation.
- Bayesian orbit determination with NUTS, using a non-centered parameterization and an optional covariance-whitened parameterization for the ill-conditioned posteriors that sparse tracking produces. Observation-noise scales are sampled, not fixed, when `infer_noise = true` (the default for the SBC study).
- Levenberg-Marquardt batch least squares with a positive-semidefinite Gauss-Newton covariance, normalized post-fit residual RMS, and a convergence flag tied to the optimizer's own stopping rule.
- Sequential EKF and UKF with a variational state-transition matrix for the EKF covariance prediction; the STM is checked against finite differences of the propagated state in `test/runtests.jl`.
- Mahalanobis D^2, marginal interval coverage, SBC ranks, and a chi-square uniformity test for those ranks.
- 2D Foster-Elrod collision probability integration.

## What is not implemented, or not yet what the name says

These are recorded as open work rather than advertised as features. See `ROADMAP.md` for the full list.

- `parse_ccsds_oem` does not parse CCSDS OEM files.
- The conjunction exporter emits a human-readable summary (`export_cdm_summary`), not a CCSDS 508.0-B-1 CDM.
- `solve_nasa_larc_uq_subproblem` summarizes posterior samples; it does not solve the NASA LaRC UQ challenge problem.
- The Foster-Elrod Pc is not validated against a reference implementation.
- SRP is not part of the Orekit cross-validation (which covers J2 and drag) and has no attitude, thermal re-radiation, or antenna-thrust model; it is not suitable for precision SRP work.
- The real-data SLR study covers one satellite (Starlette) at one station over 2.7 days, with a reduced force model; it is a validation case, not an operational POD pipeline.

## Installation

    Pkg.add(url = "https://github.com/ArpanC6/BayesianOrbitDetermination.jl")

Julia 1.10 or newer is required.

## Quick start

Two passes over one station, about 1.6 hours apart: the sparse-tracking regime this package is designed for.

    using BayesianOrbitDetermination, LinearAlgebra, Statistics

    kepler = KeplerianState(6900.0, 0.002, deg2rad(52.0), deg2rad(30.0), deg2rad(40.0), 0.0)
    cart = kepler_to_cartesian(kepler)
    u0_true = [cart.r...; cart.v...]

    gs = GroundStation("Goldstone", 35.0, -115.0, 0.5)
    windows = [(39720.0, 40110.0), (45660.0, 46050.0)]
    t_obs = vcat([collect(range(a, b, length = 8)) for (a, b) in windows]...)
    opts = OrbitPropagatorOptions(include_j2 = true)

    data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts = opts)

    guess = u0_true + [0.3, -0.2, 0.1, 3e-4, -2e-4, 1e-4]
    bls = fit_batch_least_squares(data.t_obs, data.range, data.doppler, gs,
                                  guess, 0.010, 1e-4, opts)

    C = Symmetric(Matrix(bls.covariance))
    L = Matrix(cholesky(C).L)
    chain = fit_bayesian_od(data.t_obs, data.range, data.doppler, gs,
                            bls.u_estimated, sqrt.(diag(C)), opts;
                            n_samples = 500, n_adapt = 300,
                            whitening_covariance = Matrix(C), prior_scale = 10.0)

    z = hcat([vec(chain[k]) for k in [:z1, :z2, :z3, :z4, :z5, :z6]]...)
    u_samples = (bls.u_estimated .+ 10.0 .* (L * z'))'

    ev = EvaluateCovarianceRealism(u0_true, u_samples, cov(u_samples))
    println("D^2 = ", ev.mahalanobis_d2, " (ideal 6)")
    println("90% coverage = ", ev.coverage_90)

## Tutorials

| Tutorial | Topic | Status |
|---|---|---|
| 01 | LEO J2 Bayesian orbit determination | runs; two-pass, whitened sampler |
| 02 | Drag and SRP estimation | runs; two-pass, whitened sampler |
| 03 | Batch least squares vs Bayesian NUTS | runs; includes BLS diagnostics and whitened sampler |
| 04 | TLE parsing and SGP4 propagation | runs; real SGP4, checked against a regression vector |
| 05 | EKF and UKF sequential filtering | runs; EKF covariance uses the variational STM |
| 06 | Conjunction Pc and summary export | runs; export is explicitly a summary, not a CDM |
| 07 | NASA LaRC UQ challenge | placeholder, not the challenge problem |
| 08 | Covariance realism study | runs the benchmark below |

## Tests

    julia --project=. -e 'using Pkg; Pkg.test()'

The suite covers the propagator, the observation model and its epoch convention, batch least squares, the whitened Bayesian sampler, covariance realism on synthetic Gaussian draws, the SBC uniformity test, and an end-to-end SBC simulation.

## Citation

    @software{Chakraborty_BayesianOrbitDetermination_2026,
      author = {Chakraborty, Arpan},
      title = {BayesianOrbitDetermination.jl},
      url = {https://github.com/ArpanC6/BayesianOrbitDetermination.jl},
      year = {2026}
    }

## License

MIT