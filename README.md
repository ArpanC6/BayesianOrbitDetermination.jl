# BayesianOrbitDetermination.jl

Open-Source Bayesian Orbit Determination and Covariance Realism Suite for Julia and SciML.

[![Build Status](https://github.com/ArpanC6/BayesianOrbitDetermination.jl/workflows/CI/badge.svg)](https://github.com/ArpanC6/BayesianOrbitDetermination.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

BayesianOrbitDetermination.jl is a Julia package for probabilistic orbit determination and covariance realism validation built on DifferentialEquations.jl and Turing.jl.

Traditional flight dynamics engines such as Orekit and MONTE rely on frequentist linearized Batch Least Squares, which produces overconfident state covariances under sparse ground station tracking. This package provides full Bayesian posterior sampling via NUTS MCMC with non-centered parameterization, recovering statistically calibrated uncertainty compliant with NASA NTRS UQ validation standards.

---

## Key Result: Covariance Realism Benchmark

A 10-replicate Monte Carlo benchmark on LEO J2 orbits with 50 range/Doppler observations per pass:

| Metric | Frequentist BLS | Bayesian MCMC | Ideal |
|---|---|---|---|
| Mahalanobis D^2 | 897.52 | 4.70 | 6.00 |
| 90% CI Coverage | 70% | 100% | 90% |
| Position MAE | 7.67 km | 0.28 km | 0.00 km |

The Frequentist BLS covariance is 150x overconfident (D^2 = 898 versus ideal 6), yielding only 70% interval coverage where 90% is required. The Bayesian MCMC posterior recovers statistically calibrated uncertainty (D^2 = 4.7, 100% coverage).

---

## Conjunction Assessment Impact

Covariance realism directlyC directly affects collision probability calculations. In a conjunction event at 150 m miss distance:

| Covariance Type | Position Uncertainty | Collision Probability | Compared to Realistic |
|---|---|---|---|
| Overconfident BLS | 50 m | 4.09e-3 | 1.9x overestimate |
| Realistic Bayesian | 200 m | 2.16e-3 | Calibrated |
| Bloated Conservative | 500 m | 3.91e-4 | 5x underestimate |

Overconfident covariance overestimates collision risk, triggering unnecessary avoidance maneuvers costing millions in fuel. Bloated covariance underestimates risk, potentially missing real collision threats. Only calibrated Bayesian covariance yields correct decisions.

---

## Features

- Orbital Dynamics: Two-body, J2 oblateness, J3-J5 zonal harmonics, atmospheric drag, SRP
- Observation Models: Range, Range-Rate Doppler, Azimuth, Elevation with ECEF/ECI conversion
- Bayesian Inference: NUTS MCMC via Turing.jl with non-centered parameterization for geometric efficiency
- Frequentist Baseline: Batch Least Squares (L-BFGS) with Hessian-based covariance
- Sequential Filters: Extended Kalman Filter and Unscented Kalman Filter
- Covariance Realism: Mahalanobis D^2 chi-square GOF test, 90% CI coverage, Simulation-Based Calibration
- Conjunction Assessment: 2D Foster-Elrod collision probability with CCSDS CDM export
- SGP4/TLE: TLE parser and secular J2 analytical propagator
- NASA LaRC UQ: Mixed epistemic/aleatory uncertainty challenge problem solver

---

## Installation

    Pkg.add(url = "https://github.com/ArpanC6/BayesianOrbitDetermination.jl")

---

## Quick Start

    using BayesianOrbitDetermination

    kepler = KeplerianState(6878.137, 0.001, deg2rad(51.6), deg2rad(30.0), deg2rad(40.0), deg2rad(0.0))
    cart = kepler_to_cartesian(kepler)
    u0_true = [cart.r...; cart.v...]

    gs = GroundStation("Station", 45.0, 57.5, 0.1)
    t_obs = collect(range(0.0, stop=5400.0, length=50))
    opts = OrbitPropagatorOptions(include_j2=true)

    data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts=opts)

    prior_mean = u0_true + [0.05, -0.05, 0.03, 0.0001, -0.0001, 0.0001]
    prior_std  = [0.05, 0.05, 0.05, 0.0005, 0.0005, 0.0005]

    chain = fit_bayesian_od(data.t_obs, data.range, data.doppler, gs, prior_mean, prior_std, opts)

    z_samples = hcat([vec(chain[k]) for k in [:z1, :z2, :z3, :z4, :z5, :z6]]...)
    u_samples = z_samples .* prior_std' .+ prior_mean'

    eval = EvaluateCovarianceRealism(u0_true, u_samples, cov(u_samples))
    println("Mahalanobis D^2: ", eval.mahalanobis_d2, " (Ideal = 6.0)")
    println("90% Coverage: ", eval.coverage_90)

---

## Tutorials

| Tutorial | Topic |
|---|---|
| 01 | LEO J2 Bayesian Orbit Determination |
| 02 | Drag + SRP Orbit Estimation |
| 03 | Frequentist BLS vs Bayesian MCMC Comparison |
| 04 | TLE Parsing and SGP4 Propagation |
| 05 | EKF vs UKF Sequential Filtering |
| 06 | Conjunction Risk Pc and CCSDS CDM Export |
| 07 | NASA LaRC UQ Challenge Problem |

---

## Citation

    @software{Chakraborty_BayesianOrbitDetermination_2026,
      author = {Chakraborty, Arpan},
      title = {BayesianOrbitDetermination.jl},
      url = {https://github.com/ArpanC6/BayesianOrbitDetermination.jl},
      year = {2026}
    }

## License

MIT