# BayesianOrbitDetermination.jl

Open-Source Bayesian Orbit Determination and Covariance Realism Suite for Julia and SciML.

[![Build Status](https://github.com/ArpanC6/BayesianOrbitDetermination.jl/workflows/CI/badge.svg)](https://github.com/ArpanC6/BayesianOrbitDetermination.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

`BayesianOrbitDetermination.jl` is a high-performance Julia package built on top of the SciML stack (`DifferentialEquations.jl`) and `Turing.jl` for probabilistic orbit determination, atmospheric drag / solar radiation pressure (SRP) parameter estimation, and covariance realism validation.

While traditional flight dynamics engines (such as Orekit or MONTE) rely on frequentist linearized Batch Least Squares (which frequently produce overconfident covariance estimates under sparse observation), `BayesianOrbitDetermination.jl` provides full posterior sampling with statistically calibrated covariance realism compliant with NASA NTRS UQ validation standards.

---

## Key Features

- **Orbital Dynamics**: Two-body central dynamics, $J_2$ oblateness, exponential atmospheric drag, and Solar Radiation Pressure (SRP).
- **Ground Station Observation Models**: ECEF/ECI coordinate frame conversions, Range, Range Rate (Doppler), Azimuth, and Elevation.
- **Probabilistic State Estimation**: Full Bayesian posterior sampling via No-U-Turn Sampler (NUTS) MCMC.
- **Covariance Realism Suite**: Mahalanobis distance $\chi^2$ Goodness-of-Fit (GOF) testing, Empirical Cumulative Distribution Function (ECDF) calibration, and Simulation-Based Calibration (SBC).
- **Baseline Comparison**: Built-in frequentist Batch Least Squares (Gauss-Newton / L-BFGS) baseline.

---

## Preliminary Benchmark Summary (100 Replicates)

Across a 100-trial Monte Carlo benchmark on Low Earth Orbit (LEO) trajectories:

| Metric | Frequentist Batch Least Squares | Bayesian MCMC (BayesianOrbitDetermination.jl) | Ideal |
|---|---|---|---|
| **Mahalanobis Distance $D^2$** | $14.82 \pm 0.95$ | $6.24 \pm 0.31$ | $6.00$ |
| **Covariance Realism (GOF Pass)** | $32\%$ | $94\%$ | $95\%$ |
| **90% Interval Coverage** | $68.5\%$ | $91.2\%$ | $90.0\%$ |

Full methodology and per-parameter standard error breakdowns are documented in `benchmarks/RESULTS.md`.

---

## Installation

`BayesianOrbitDetermination.jl` can be installed directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/ArpanC6/BayesianOrbitDetermination.jl")
```

---

## Quick Start Example

```julia
using BayesianOrbitDetermination
using LinearAlgebra

# 1. Define LEO satellite state & ground station
kepler = KeplerianState(6878.137, 0.001, deg2rad(51.6), deg2rad(30.0), deg2rad(0.0), deg2rad(0.0))
cart = kepler_to_cartesian(kepler)
u0_true = [cart.r...; cart.v...]

gs = GroundStation("Hartebeesthoek", -25.887, 27.707, 1.558)
t_obs = collect(range(0.0, stop=3600.0, length=20))

# 2. Simulate synthetic tracking observations (Range & Doppler)
opts = OrbitPropagatorOptions(include_j2=true)
data = simulate_tracking_data(u0_true, t_obs, gs, 0.010, 1e-4; opts=opts)

# 3. Fit Bayesian Orbit Determination
prior_mean = u0_true + [1.0, -1.0, 0.5, 0.001, -0.001, 0.001]
prior_std  = [5.0, 5.0, 5.0, 0.005, 0.005, 0.005]

chain = fit_bayesian_od(t_obs, data.range, data.doppler, gs, prior_mean, prior_std, opts)

# 4. Assess Covariance Realism
u_samples = Matrix(chain[:u0])
cov_est = cov(u_samples)
realism = EvaluateCovarianceRealism(u0_true, u_samples, cov_est)

println("Mahalanobis Distance D^2: ", realism.mahalanobis_d2)
println("Chi2 GOF p-value: ", realism.p_value_gof)
```

---

## Repository Structure

```
BayesianOrbitDetermination.jl/
├── Project.toml
├── README.md
├── ROADMAP.md
├── LICENSE
├── CITATION.cff
├── CONTRIBUTING.md
├── src/
│   ├── BayesianOrbitDetermination.jl
│   ├── dynamics/
│   │   ├── two_body.jl
│   │   ├── perturbations.jl
│   │   └── propagator.jl
│   ├── observations/
│   │   └── measurement_models.jl
│   ├── inference/
│   │   ├── bayesian_od.jl
│   │   └── optimizers.jl
│   └── realism/
│       ├── covariance_realism.jl
│       └── sbc.jl
├── examples/
│   ├── 01_two_body_j2_od.jl
│   ├── 02_drag_srp_estimation.jl
│   └── 03_frequentist_vs_bayesian_od.jl
├── test/
│   └── runtests.jl
├── benchmarks/
│   ├── run_covariance_realism_study.jl
│   └── RESULTS.md
└── .github/
    └── workflows/
        └── CI.yml
```

---

## Citation

If you use `BayesianOrbitDetermination.jl` in academic work or agency research, please cite:

```bibtex
@software{Chakraborty_BayesianOrbitDetermination_2026,
  author = {Chakraborty, Arpan},
  title = {BayesianOrbitDetermination.jl: Open-Source Bayesian Orbit Determination and Covariance Realism Suite},
  url = {https://github.com/ArpanC6/BayesianOrbitDetermination.jl},
  year = {2026}
}
```

---

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
