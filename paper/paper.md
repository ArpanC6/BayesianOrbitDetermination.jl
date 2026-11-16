---
title: 'BayesianOrbitDetermination.jl: Bayesian Orbit Determination and Covariance Realism for Sparse Tracking'
tags:
  - Julia
  - Astrodynamics
  - Orbit Determination
  - Bayesian Inference
  - Uncertainty Quantification
  - Covariance Realism
authors:
  - name: Arpan Chakraborty
    orcid: 0009-0006-7735-8520
    affiliation: 1
affiliations:
  - name: Department of Computer Science and Engineering, Sikkim Global Technical University
    index: 1
date: 16 November 2026
bibliography: paper.bib
---

# Summary

`BayesianOrbitDetermination.jl` is an open-source Julia package for probabilistic orbit determination (OD) and covariance realism validation. It couples `DifferentialEquations.jl` [@rackauckas2017differentialequations] with the No-U-Turn Sampler (NUTS) in `Turing.jl` [@ge2018turing] to produce full posterior distributions over spacecraft states, alongside a converged Levenberg–Marquardt batch least-squares (BLS) estimator with Gauss–Newton covariance and an explicit convergence flag.

The package addresses a specific gap in open-source astrodynamics: existing tools report point estimates with linearized covariances, but provide no built-in diagnostic to test whether those covariances are calibrated under sparse tracking conditions. `BayesianOrbitDetermination.jl` provides Mahalanobis $\chi^2$ goodness-of-fit tests, marginal interval coverage, simulation-based calibration (SBC) ranks with uniformity testing [@talts2018validating], and a covariance-whitened sampling option for ill-conditioned posteriors ($\kappa \sim 10^{10}$).

# Statement of need

Space Traffic Management (STM) and conjunction assessment depend on realistic covariance matrices [@poore2016covariance]. Overconfident covariances lead to missed collision warnings, inflated covariances cause unnecessary maneuvers. Operational flight dynamics suites such as Orekit [@orekit] and MONTE rely primarily on frequentist batch least squares with linearized covariances, but do not include tools to test whether those covariances are statistically calibrated.

Under sparse ground-station tracking - a few short passes separated by hours - the state covariance becomes severely ill-conditioned. Standard samplers stall in state coordinates, and frequentist optimizers may mislabel convergence. The result is that comparisons between linearized and Bayesian covariances are often confounded by computation failures rather than reflecting a property of the methods.

`BayesianOrbitDetermination.jl` is built so that both estimators can be computed correctly on identical data: a converged BLS fit and a whitened NUTS posterior. The target audience is researchers and practitioners in astrodynamics and space situational awareness who need to evaluate covariance realism for sparse-tracking scenarios, and who want reproducible, open-source tools for that evaluation.

# State of the field

Open-source astrodynamics packages with orbit determination capabilities include Orekit (Java) [@orekit], `Astrodynamics.jl` (Julia) [citation:5], and `EphemerisSources.jl` (Julia) [citation:10]. Orekit is the most widely used open-source flight dynamics library and provides BLS estimation with linearized covariance. It does not provide Bayesian posterior sampling or covariance calibration diagnostics.

`Astrodynamics.jl` provides a Julia framework for orbital mechanics, including propagation and optimization. It focuses on trajectory design and optimization rather than statistical inference on OD posteriors. `EphemerisSources.jl` addresses ephemeris sourcing and parsing, not orbit determination.

A separate contribution of this package is the whitened parameterization for ill-conditioned posteriors. Sampling in state coordinates fails when the covariance condition number reaches $10^{10}$, which is typical for sparse LEO tracking with a single station. The whitening transform converts the stiff, anisotropic posterior into an isotropic space, enabling NUTS to achieve zero divergences and effective sample sizes of order the number of draws. This design decision is described in the Software design section.

# Software design

**Covariance whitening.** The central design decision is the parameterization for NUTS. Rather than sampling directly in state coordinates, the package samples in a whitened space:

$$\mathbf{u}_0 = \mathbf{u}_{\text{BLS}} + s\,L\mathbf{z}, \qquad LL^T = P_{\text{BLS}}, \qquad \mathbf{z} \sim \mathcal{N}(0, I_6),$$

where $\mathbf{u}_{\text{BLS}}$ and $P_{\text{BLS}}$ come from a converged BLS fit, and $s$ controls prior width in units of linearized standard deviations. This transform is the difference between a sampler that stalls and one that achieves zero divergences. The trade-off is that the prior is centered on the BLS estimate, making the posterior empirical-Bayes rather than prior-blind.

**Dual estimators on identical data.** BLS and Bayesian estimation share the same force model, measurement model, and data. This ensures any comparison reflects the inference method, not differences in modeling. The BLS implementation includes a convergence flag tied to the optimizer’s stopping rule, without this, an unconverged fit can make the linearized covariance appear badly overconfident, which is an optimizer failure rather than a property of the Gauss–Newton approximation.

**Measurement model and epoch handling.** For Satellite Laser Ranging (SLR), the CRD format reports ground-transmit epochs (epoch event 2). The parser converts each epoch to reception time before the two-way light-time solution. Applying the range model without this conversion inflates the CPF residual from 7.83 m to 39.5 m RMS. This convention was identified by diagnostic, not assumed.

**What is not implemented.** The conjunction exporter emits a human-readable summary, not a CCSDS 508.0-B-1 CDM. `parse_ccsds_oem` does not parse CCSDS OEM files. The Foster - Elrod Pc is not validated against a reference implementation. SRP is not part of the Orekit cross-validation. These gaps are recorded in `ROADMAP.md` rather than hidden behind feature names.

# Research impact statement

**Independent validation on real SLR data.** The package ingests 98 ILRS CRD v2 normal points to Starlette from Yarragadee (station 7090), 2025-06-09 to 2025-06-11, and agrees with the independent ILRS CPF ephemeris to **7.83 m RMS** (mean 4.15 m, max 19.1 m) over 2.74 days at 15–66° elevation [@pearlman2019ilrs; @ilrs_crd]. The light-time correction is 20.9 m RMS. This is a measurement-model validation on real tracking data, not a self-consistent simulation.

**Industry-standard cross-validation.** The numerical propagator is cross-validated against CNES Orekit 13.1.9: maximum position difference over 24 h is **1.2 mm for J2** and **2.8 mm for J2 + drag**. This establishes agreement with the most widely used open-source flight dynamics library at millimetre level for the force models covered.

**Covariance realism benchmark.** A 50-scenario study produces `benchmarks/RESULTS.md` with per-trial Mahalanobis $D^2$, coverage, convergence flags, divergences, $\hat{R}$, and ESS. The mean $D^2$ is $7.37 \pm 0.74$ for BLS and $7.89 \pm 0.80$ for the posterior (ideal: 6). Per-axis position MAE is statistically indistinguishable (26.3 m vs 26.2 m). The SBC study (N=104) reports zero divergences and worst $\hat{R} = 1.046$. These are reproducible benchmark materials generated by committed scripts, not hand-typed values.

# AI usage disclosure

Generative AI tools were used during the preparation of this manuscript and the associated software. Specifically, AI was used for code refactoring suggestions, documentation drafting, and editorial review of the manuscript text. All software design decisions, implementation, validation, and scientific interpretation were performed and verified by the author. All AI-generated suggestions were reviewed, modified, and tested by the author before inclusion. The author takes full responsibility for the accuracy and correctness of all submitted materials.

# Acknowledgements

The author acknowledges the SciML open-source community and the International Laser Ranging Service (ILRS) for making tracking data openly accessible.

# References
