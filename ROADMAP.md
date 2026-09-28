# Development Roadmap

## Project Vision

To become the open-source gold standard for statistically calibrated Bayesian orbit determination, covariance realism verification, and uncertainty quantification across academic research, aerospace industry, and space agency applications including NASA, ESA, and CNES.

---

## Phase 1: Core Engine and Covariance Realism (Months 1 to 3) - COMPLETED

- [x] Keplerian and Cartesian Dynamics Engine: ECI state transformations and two-body propagation
- [x] Orbital Perturbations: J2 zonal harmonic, exponential atmospheric drag, Solar Radiation Pressure
- [x] Ground Station Observation Models: ECEF-to-ECI rotations, Range, Range-Rate Doppler, Azimuth, Elevation
- [x] Bayesian Inference Engine: Turing.jl model with NUTS MCMC and non-centered parameterization
- [x] Covariance Realism Suite: Mahalanobis D^2 chi-square GOF test, 90% CI coverage, Simulation-Based Calibration
- [x] 10-Replicate Benchmark Study: Monte Carlo validation comparing Frequentist BLS against Bayesian MCMC

Key result: Frequentist BLS D^2 = 897.52 (150x overconfident), Bayesian MCMC D^2 = 4.70 (calibrated)

---

## Phase 2: TLE/SGP4 Integration and Sequential Filtering (Months 4 to 6) - COMPLETED

- [x] SGP4 Analytical Propagator Bridge: TLE Line 1/Line 2 parser and secular J2 drift propagator
- [x] CCSDS OEM Tracking Data Ingestion: Ground Station Observation parsing module
- [x] High-Order Gravitational Spherical Harmonics: J2, J3, J4, J5 zonal harmonic acceleration model
- [x] Sequential Kalman Filtering: EKF and UKF with ForwardDiff autodiff Jacobians and matrix exponential STM

---

## Phase 3: Conjunction Assessment and Collision Risk (Months 7 to 9) - COMPLETED

- [x] Probability of Collision Engine: 2D Foster-Elrod Pc integration over encounter plane
- [x] Space Traffic Management CCSDS Export: Autonomous Conjunction Data Message generator
- [x] Overconfidence Impact Assessment: Conjunction risk case study showing covariance choice affects Pc by 10x

Key result: Overconfident BLS overestimates Pc by 1.9x, bloated covariance underestimates by 5x

---

## Phase 4: Agency Benchmarking, JOSS Paper and Python Interop (Months 10 to 12) - COMPLETED

- [x] NASA LaRC UQ Challenge Solver: Mixed epistemic/aleatory uncertainty bound solver
- [x] JOSS Manuscript: Peer-review publication manuscript (paper/paper.md and paper/paper.bib)
- [x] Python Interop Interface: Python wrapper using juliacall for Orekit, Poliastro, Astropy integration

---

## Future Work

- [ ] 100-replicate full benchmark study with standard error analysis
- [ ] Higher-order gravitational models (J6+) and third-body perturbations
- [ ] Multi-object conjunction screening pipeline
- [ ] Integration with Orekit for validated dynamics comparison
- [ ] GPU-accelerated MCMC via CUDA.jl for large constellation OD
- [ ] Extended SGP4 with deep-space corrections