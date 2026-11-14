# Development Roadmap

## Current state

Every item below was checked against the code, not against intent.

| Area | State |
|---|---|
| Two-body + J2 + drag + SRP propagation | Implemented; exact inverse-square gravity; numerical two-body matches the analytic Kepler solution to millimetres over 24 h at tight tolerances; J2 and J2+drag match Orekit 13.1.9 to 1.2 mm and 2.8 mm over 24 h; J2 nodal drift matches the analytic secular rate within 1% over 50 orbits; SRP uses an analytic solar ephemeris and a conical shadow, but is not part of the Orekit cross-validation |
| Zonal harmonics | Implemented for J2-J5 by Legendre recurrence, wired into the propagator behind `include_j3j4 = true`; checked against the J2 closed form and against a finite-difference gradient of the potential in the test suite |
| Range / Doppler observation model | Implemented; WGS-84 geodetic to ECEF to ECI; two-way light time solved by fixed-point iteration (dominant correction for SLR, ~20 m at LEO ranges) |
| Bayesian OD (NUTS) | Implemented; non-centered and covariance-whitened parameterizations; noise scales inferred from the data with `infer_noise = true` |
| Batch least squares | Implemented with Levenberg-Marquardt; Gauss-Newton covariance; returns residual RMS, iteration count, Newton decrement and a convergence flag tied to the optimizer's stopping rule; range-only mode for SLR |
| Covariance realism diagnostics | Mahalanobis D^2, marginal coverage, SBC ranks with a uniformity test and rank-histogram figure, and a sharded SBC study |
| Foster-Elrod Pc | Implemented; not validated against a reference implementation |
| TLE parsing and SGP4 | Real SGP4 (SatelliteToolboxSgp4, Vallado); checked against a published regression vector, not Kepler + J2 |
| CCSDS OEM ingestion | Not implemented; current parser reads a different format |
| CCSDS CDM export | Not implemented; the exporter is named `export_cdm_summary` and labels itself as a summary |
| EKF / UKF | Implemented; EKF covariance prediction uses the variational STM, checked against finite differences. Not validated against flight data |
| NASA LaRC UQ challenge | Not implemented |

## Next: make the existing claims true

1. ~~Propagator validation against an independent reference (Orekit, same
   scenario, same force model) at the 1 m level over one day.~~ Done for J2 and
   J2+drag: Orekit 13.1.9 agrees to 1.2 mm and 2.8 mm over 24 h
   (`validation/orekit/`). SRP is not part of that comparison.
2. ~~Real SGP4, or rename `propagate_sgp4` to what it does.~~ Done: SGP4 comes
   from SatelliteToolboxSgp4 (Vallado) and is tested against a regression vector.
3. ~~Wire J3/J4 into `orbit_dynamics` behind a flag; add J5 or drop the claim.~~
   Done: `include_j3j4 = true` adds J2-J5 from the Legendre recurrence.
4. ~~Implement the EKF state-transition matrix through the variational
   equations.~~ Done; `propagate_state_and_stm` integrates `dPhi/dt = (df/du) Phi`
   with the Jacobian from forward-mode AD, and a test compares it against finite
   differences of the propagated state. Cross-validating the EKF against flight
   data or an independent filter is still open.
5. ~~Implement CCSDS 508.0-B-1 CDM output, or rename the exporter.~~ Renamed to
   `export_cdm_summary`; the file says it is not a CDM. Conformant output is
   still open and needs a calendar TCA plus 6x6 RTN covariances.
6. Implement the NASA LaRC challenge subproblems, or remove the module.

## Next: make the statistics honest

7. ~~Infer observation-noise scales in the Bayesian model instead of hard-coding them.~~
   Done: `sr` and `sd` are sampled in `bayesian_orbit_determination_model`
   whenever `infer_noise = true`, and the SBC study treats them as parameters.
8. Report NUTS diagnostics (r-hat, ESS, divergences) for every chain that enters
   a result. (Done in the covariance study; `results_per_trial.csv` records them
   per replication.)
9. ~~Turn the SBC ranks into a uniformity test with a figure.~~ Done; see
   `benchmarks/run_sbc_study.jl`, `benchmarks/SBC.md`, and
   `benchmarks/sbc_rank_histogram.svg`.
10. Regenerate `benchmarks/RESULTS.md` from a single script, with seeds and
    convergence filters recorded alongside the numbers. (Done; the script is
    `benchmarks/run_covariance_realism_study.jl`.)

## What the covariance study has established so far

`benchmarks/RESULTS.md` is the authoritative record; this is the interpretation.

- With a converged Levenberg-Marquardt fit, the Gauss-Newton covariance from
  `(J'J)^-1` is approximately calibrated on self-consistent simulated data. In
  the definitive 50-trial run the mean Mahalanobis D^2 is 7.37 +/- 0.74, 1.9
  standard errors above 6, and the per-component 90% intervals cover the truth
  at 90-94%. The Bayesian mean is 7.89 +/- 0.80 and the paired Bayes-minus-BLS
  difference is +0.52 +/- 0.20, so the full posterior inherits the residual
  rather than removing it. Whether the residual is a small real overconfidence
  in the least-constrained direction or seed variance needs more replications
  than the 50 recorded here.
- The earlier claim that the linearized covariance is "150x overconfident" was
  an artifact of an optimizer that had not converged, not a property of the
  Gauss-Newton covariance.
- Sparse tracking makes the posterior covariance ill-conditioned. In the sparse
  scenarios examined here the condition number is of order 1e10, and the
  directions with the smallest eigenvalues are combinations of along-track and
  radial components. A sampler in state coordinates cannot move along those
  directions with a sane step size.
- Whitening by the BLS covariance, `u = u_BLS + s L z` with `L L' = P_BLS`,
  makes the geometry isotropic and lets NUTS sample the same posterior with
  order-unit effective sample sizes. This is the parameterization the benchmark
  and Tutorial 03 use.
- On this problem, once both estimators are computed correctly, the Bayesian
  posterior and the linearized covariance agree on the state and are both
  approximately calibrated: marginal coverage is nominal, while the joint D^2
  carries the mild positive bias described above for both estimators. The
  interesting Bayesian advantage, if it exists, is not in the point estimate or
  the marginal covariance. It has to be looked for in non-Gaussianity of the
  posterior, in model error, or in cases where the linearization is genuinely
  invalid.
- A strict absolute Newton-decrement test mislabels converged fits. In 15 random
  sparse-arc scenarios, restarting the optimizer with 500 iterations moved the
  estimate by exactly zero in all 15, but the absolute test flagged 1 of them as
  unconverged. The cause is dividing a rounding-level gradient by the smallest
  eigenvalue of `J'J`, which is ~1e-10 for these arcs. The convergence flag now
  follows the optimizer's own stopping rule, and `newton_decrement` is returned
  as a diagnostic instead of deciding the flag.

## Open research question

Whether the commonly asserted advantage of full Bayesian posteriors over
linearized batch least squares — better-calibrated covariance under sparse
tracking — survives when the optimizer is actually converged, the noise scales
are inferred rather than fixed, and calibration is assessed by SBC instead of a
single Mahalanobis statistic.

The definitive study answers the fixed-noise, D^2-and-coverage version of the
question with a negative: no detectable advantage, and the posterior is
$+0.52\pm0.20$ *worse* in mean D^2 than the BLS covariance. The inferred-noise
and SBC versions remain open, and they are the experiments this package should
be used for next.
