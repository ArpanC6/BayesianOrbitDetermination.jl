---
title: 'BayesianOrbitDetermination.jl: Open-Source Bayesian Orbit Determination and Covariance Realism Suite'
tags:
  - Julia
  - Astrodynamics
  - Orbit Determination
  - Scientific Machine Learning
  - Uncertainty Quantification
  - Covariance Realism
authors:
  - name: Arpan Chakraborty
    affiliation: 1
affiliations:
  - name: Department of Computer Science and Engineering, Sikkim Global Technical University
    index: 1
date: 25 October 2026
bibliography: paper.bib
---

# Summary

`BayesianOrbitDetermination.jl` is an open-source Julia package for probabilistic
orbit determination and covariance realism validation within the Scientific
Machine Learning (SciML) ecosystem. Flight dynamics engines such as Orekit
[@orekit] and MONTE rely primarily on frequentist batch least squares with
linearized covariances; whether full posterior sampling materially improves the
calibration of those covariances under sparse ground-station tracking is an open
question that this package is built to test.

The package couples numerical differential equation solvers from
`DifferentialEquations.jl` [@rackauckas2017differentialequations] with No-U-Turn
Sampler (NUTS) Hamiltonian Monte Carlo in `Turing.jl` [@ge2018turing], and
provides Mahalanobis $\chi^2$ goodness-of-fit tests, marginal interval coverage,
simulation-based calibration ranks with a chi-square uniformity test
[@talts2018validating], 2D Foster-Elrod collision probability integration
[@foster1992matrix], and sequential EKF/UKF filters whose EKF covariance
prediction integrates the variational equations rather than assuming an
identity state-transition matrix. The range model is validated against 98
International Laser Ranging Service normal points to Starlette, where the
package agrees with an independent ILRS CPF prediction to 7.83 m RMS over
2.74 days [@pearlman2019ilrs].

# Statement of Need

Space Traffic Management (STM) and satellite collision risk assessment depend on
realistic covariance matrices [@poore2016covariance]. Overconfident covariance
estimates lead to missed collision warnings, whereas inflated covariances cause
excessive operator maneuvers. `BayesianOrbitDetermination.jl` provides full
posterior distributions together with the diagnostics needed to test whether they
are calibrated, and a conjunction-data export that is clearly distinguished from
CCSDS 508.0-B-1 conformance.

# Mathematics and Dynamics Models

The state transition of a satellite state $\mathbf{x}(t) = [\mathbf{r}(t)^T, \mathbf{v}(t)^T]^T$ is governed by central gravitational attraction and perturbations [@vallado2013fundamentals; @tapley2004statistical]:

$$\ddot{\mathbf{r}} = -\frac{\mu}{\|\mathbf{r}\|^3}\mathbf{r} - \nabla\!\left[\frac{\mu}{r}\sum_{n=2}^{5} J_n \left(\frac{R_\oplus}{r}\right)^n P_n\!\left(\frac{z}{r}\right)\right] + \mathbf{a}_{\text{drag}} + \mathbf{a}_{\text{SRP}},$$

where the zonal sum is evaluated with a Legendre recurrence. The observational measurement likelihood for range $\rho$ and range rate $\dot{\rho}$ observed from a ground station is modeled as:

$$y_{\text{range}} \sim \mathcal{N}\left(\|\mathbf{r} - \mathbf{r}_{\text{gs}}\|, \sigma_{\text{range}}^2\right), \quad y_{\text{doppler}} \sim \mathcal{N}\left(\frac{(\mathbf{r} - \mathbf{r}_{\text{gs}}) \cdot (\mathbf{v} - \mathbf{v}_{\text{gs}})}{\|\mathbf{r} - \mathbf{r}_{\text{gs}}\|} , \sigma_{\text{doppler}}^2\right)$$

For laser ranges the observation is the two-way travel time, and the reported
range refers to the satellite at the bounce epoch. The model therefore
evaluates the satellite at that epoch and rotates the station through the up-
and down-leg light times, solving the coupled travel times by fixed-point
iteration. At Starlette's roughly 1300 km range the bounce epoch is 4.3 ms
before reception and the correction changes the range by 20.9 m RMS.

# Sampling an Ill-Conditioned OD Posterior

Sparse tracking, a few short passes separated by hours, produces a posterior
whose covariance has a condition number of order $10^{10}$ in the scenarios
examined here. The stiff directions mix along-track position with range-rate
information, and a sampler in state coordinates cannot take a step that is
simultaneously sane in the well- and ill-determined directions. The package
therefore offers a covariance-whitened parameterization,

$$\mathbf{u}_0 = \mathbf{u}_{\text{BLS}} + s\,L\mathbf{z}, \qquad LL^T = P_{\text{BLS}}, \qquad \mathbf{z} \sim \mathcal{N}(0, I_6),$$

where $\mathbf{u}_{\text{BLS}}$ and $P_{\text{BLS}}$ come from a converged batch
least-squares fit and $s$ sets the prior width in units of the linearized
standard deviations. Sampling in $\mathbf{z}$ is isotropic and lets NUTS reach
effective sample sizes of the order of the number of draws in the benchmark
scenarios.

# Validation and Reproducibility

No result is typed into the documentation by hand. The covariance-realism study
in `benchmarks/run_covariance_realism_study.jl` simulates 50 independent LEO
scenarios with sparse range and Doppler tracking, fits the same data by
Levenberg-Marquardt batch least squares and by whitened NUTS, and writes
`benchmarks/RESULTS.md` together with per-trial values, convergence flags,
divergences, $\hat{R}$, and effective sample sizes in
`benchmarks/results_per_trial.csv`.

Two properties are checked: the Mahalanobis $D^2$ of the true state against the
reported covariance, whose expectation is 6 for a six-dimensional state, and the
fraction of trials in which the truth falls inside each marginal 90% interval.
Under this protocol the two estimators produce statistically indistinguishable
position errors (26.3 m and 26.2 m per-axis MAE), and the per-component 90%
intervals cover the truth at 90-94%. The mean
Mahalanobis $D^2$ is $7.37\pm0.74$ for the linearized covariance and
$7.89\pm0.80$ for the posterior, against an ideal value of 6; the paired
difference is $+0.52\pm0.20$. The full posterior therefore does not remove the
residual optimism, and 50 replications are not enough to decide whether it is a
small real bias or seed variance. The study also records the negative result
that motivated the package's whitening option: an unconverged fit can make the
linearized covariance look badly overconfident, which is an optimizer failure
rather than a property of the Gauss-Newton approximation.

# Independent Validation on Satellite Laser Ranging

The covariance-realism study is self-consistent: the same force and measurement
models generate and fit the data. The `validation/slr/` case adds a real-data
check using 98 ILRS CRD v2 normal points to Starlette from Yarragadee (station
7090) on 2025-06-09 to 2025-06-11 and the matching ILRS Consolidated Prediction
Format ephemeris [@pearlman2019ilrs; @ilrs_crd].

The CRD record 11 epochs in this data set are ground-transmit tagged (epoch
event 2), so each epoch is converted to reception time before the two-way
light-time solution. Applying the range model to the independent CPF prediction
gives a 7.83 m RMS residual over the 98 normal points (mean 4.15 m, maximum
19.1 m) at elevations of 15-66 degrees and a 2.74-day span. The light-time
correction is 20.9 m RMS; applying it without the epoch conversion inflates the
residual to 39.5 m, which is how the convention was identified. The remaining
+4 m mean is consistent with unmodelled tropospheric refraction.

A range-only Levenberg-Marquardt fit with the package's force model reaches
24.7 m RMS over the first day (46 normal points) and 58.8 m RMS over the full
arc. The growth is consistent with the reduced force model (no tesseral
harmonics, solid Earth tides, or atmospheric delay) and with the CPF being an
independent prediction rather than truth. This case validates the measurement
model and exercises the estimator on real data; it is not an operational POD
accuracy claim.

# Limitations

The covariance study is deliberately self-consistent: truth and estimator use the same
force model, so it tests the estimator and the diagnostics, not the force model;
the SLR case above removes that limitation for the measurement model and adds a
real-data check of the dynamics, with the force-model caveats stated there.
The validation study fixes observation-noise scales at their generating values;
the separate SBC study samples them. The Bayesian prior is centered on the batch
least-squares estimate, so the reported posterior is empirical-Bayes, not a
prior-blind posterior. The scenarios are LEO with a single ground station. The
two-body core matches the closed-form Kepler solution to millimetres over one
day, the J2 term reproduces the analytic secular nodal rate within 1% over 50
orbits, and the J2 and J2+drag propagator is cross-validated against Orekit
13.1.9 to 1.2 mm and 2.8 mm over 24 h. The J3-J5 zonal terms are implemented by
Legendre recurrence and checked against a finite-difference gradient of the
potential. The solar radiation pressure model uses an analytic solar ephemeris
and a conical shadow, but it is not included in the Orekit validation. The
conjunction export is a documented key-value summary rather than a CCSDS
508.0-B-1 CDM [@ccsds2013cdm], and SGP4 comes from SatelliteToolboxSgp4 rather
than being reimplemented here. These gaps are tracked in `ROADMAP.md` rather
than hidden behind feature names.

# Acknowledgements

The author acknowledges the SciML open-source community and Julia Space initiatives.

# References
