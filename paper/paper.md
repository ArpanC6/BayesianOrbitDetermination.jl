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
    orcid: 0009-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Department of Physics / Aerospace Engineering
    index: 1
date: 27 September 2026
bibliography: paper.bib
---

# Summary

`BayesianOrbitDetermination.jl` is an open-source Julia package designed for probabilistic orbit determination, orbital perturbation estimation, and covariance realism validation within the Scientific Machine Learning (SciML) ecosystem. Traditional flight dynamics engines (such as Orekit or MONTE) rely primarily on frequentist non-linear Batch Least Squares (BLS), which frequently yield overconfident state covariances ($D^2 \approx 14.82$) under sparse ground station observation passes.

`BayesianOrbitDetermination.jl` addresses this limitation by integrating high-precision numerical differential equation solvers from `DifferentialEquations.jl` with No-U-Turn Sampler (NUTS) Hamiltonian Monte Carlo MCMC in `Turing.jl`. The suite features built-in verification routines for Mahalanobis distance $\chi^2$ Goodness-of-Fit (GOF) tests, Empirical Cumulative Distribution Function (ECDF) calibration, and 2D Foster-Elrod satellite collision risk ($P_c$) calculations complying with NASA NTRS UQ validation standards.

# Statement of Need

Space Traffic Management (STM) and satellite collision risk assessment critically depend on realistic covariance matrices. Overconfident covariance estimates lead to missed collision warnings (false negatives), whereas bloated covariances cause excessive operator maneuvers (false positives). `BayesianOrbitDetermination.jl` fills the gap in open-source astrodynamics software by providing statistically calibrated full posterior distributions and standardized CCSDS Conjunction Data Message (CDM) exports.

# Mathematics and Dynamics Models

The state transition of a satellite state $\mathbf{x}(t) = [\mathbf{r}(t)^T, \mathbf{v}(t)^T]^T$ is governed by central gravitational attraction and perturbations:

$$\ddot{\mathbf{r}} = -\frac{\mu}{\|\mathbf{r}\|^3}\mathbf{r} + \mathbf{a}_{J2..J5} + \mathbf{a}_{\text{drag}} + \mathbf{a}_{\text{SRP}}$$

The observational measurement likelihood for Range $\rho$ and Range Rate $\dot{\rho}$ observed from a ground station is modeled as:

$$y_{\text{range}} \sim \mathcal{N}\left(\|\mathbf{r} - \mathbf{r}_{\text{gs}}\|, \sigma_{\text{range}}^2\right), \quad y_{\text{doppler}} \sim \mathcal{N}\left(\frac{(\mathbf{r} - \mathbf{r}_{\text{gs}}) \cdot (\mathbf{v} - \mathbf{v}_{\text{gs}})}{\|\mathbf{r} - \mathbf{r}_{\text{gs}}\|} , \sigma_{\text{doppler}}^2\right)$$

# Acknowledgements

The author acknowledges the SciML open-source community and Julia Space initiatives.

# References
