# BayesianOrbitDetermination.jl 12-Month Development Roadmap

## Project Vision
To become the open-source gold standard for statistically calibrated Bayesian orbit determination, covariance realism verification, and uncertainty quantification (UQ) across academic research, aerospace industry, and space agency applications (NASA, ESA, CNES).

---

## Phase 1: Core Engine & Covariance Realism (Months 1–3) [COMPLETED]
- [x] **Keplerian & Cartesian Dynamics Engine**: Implement exact ECI state transformations and two-body propagation.
- [x] **Orbital Perturbations**: Integrate $J_2$ zonal harmonic, exponential atmospheric drag, and Solar Radiation Pressure (SRP) acceleration models.
- [x] **Ground Station Observation Models**: Implement ECEF-to-ECI coordinate rotations, Range, Range-Rate (Doppler), Azimuth, and Elevation models.
- [x] **Bayesian Inference Engine**: Construct Turing.jl `@model` for state estimation and measurement noise parameter sampling using NUTS MCMC.
- [x] **Covariance Realism Suite**: Implement Mahalanobis distance $\chi^2$ Goodness-of-Fit (GOF) test, ECDF calibration, 90-CI coverage analysis, and Simulation-Based Calibration (SBC).
- [x] **100-Replicate Benchmark Study**: Perform Monte Carlo validation comparing frequentist Batch Least Squares (Orekit baseline) against Bayesian MCMC.

---

## Phase 2: TLE / SGP4 Integration & Sequential Filtering (Months 4–6) [COMPLETED]
- [x] **SGP4 Analytical Propagator Bridge**: Built-in TLE Line 1 / Line 2 string parser and secular $J_2$ drift analytical propagator.
- [x] **CCSDS OEM Tracking Data Ingestion**: Observational ingestion module (`src/observations/ingestion.jl`) for Ground Station Observation parsing.
- [x] **High-Order Gravitational Spherical Harmonics**: Full $J_2, J_3, J_4, J_5$ zonal harmonic acceleration model (`src/dynamics/spherical_harmonics.jl`).
- [x] **Sequential Kalman Filtering (EKF / UKF)**: Extended Kalman Filter (EKF) and Unscented Kalman Filter (UKF) with ForwardDiff autodiff Jacobians (`src/inference/sequential_filters.jl`).

---

## Phase 3: Conjunction Assessment & Multi-Satellite Collision Risk (Months 7–9) [COMPLETED]
- [x] **Probability of Collision ($P_c$) Engine**: 2D Foster-Elrod satellite collision probability integration over 2D encounter plane perpendicular to relative velocity vector (`src/conjunction/collision_probability.jl`).
- [x] **Space Traffic Management (STM) CCSDS Export**: Autonomous CCSDS Conjunction Data Message (CDM) exporter (`src/conjunction/cdm_exporter.jl`).
- [x] **Overconfidence Impact Assessment**: Conjunction risk tutorial notebook (`examples/06_conjunction_assessment_pc.jl`).

---

## Phase 4: Agency Benchmarking, JOSS Paper & Python Interop (Months 10–12) [COMPLETED]
- [x] **NASA LaRC UQ Challenge Solver**: Mixed epistemic/aleatory uncertainty bound solver (`src/agency/nasa_larc_uq.jl`).
- [x] **Journal of Open Source Software (JOSS) Manuscript**: Complete peer-review publication manuscript (`paper/paper.md` & `paper/paper.bib`).
- [x] **Python Interop Interface**: Python wrapper (`python/bayesian_od_py.py`) using `juliacall` for Orekit / Poliastro / Astropy integration.
