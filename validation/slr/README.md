# Real-data validation: satellite laser ranging to Starlette

The benchmark and SBC studies are self-consistent: the same force model and
observation model generate and fit the data. This directory removes that
limitation for one case by fitting **real ILRS laser-ranging data** to a
satellite whose orbit is known from an independent ILRS product.

`run_starlette.jl` performs three checks:

1. **CPF prediction vs normal points.** Interpolate the ILRS Consolidated
   Prediction Format ephemeris at the reception epoch, compute the two-way
   range through the same measurement model, and compare with the laser
   normal points.
2. **Our force model vs normal points.** Start from the CPF state and propagate
   with our J2-J5 + drag + SRP model.
3. **Range-only batch least squares.** Estimate the six-component state at the
   first normal point and compare the post-fit residuals.

## Data

| Item | Value |
|---|---|
| Satellite | Starlette (ILRS 7501001, NORAD 07646) |
| Shape / mass | 24 cm sphere, 47.3 kg, A/m = 9.57e-4 m^2/kg |
| Station | Yarragadee, Australia (YARL, 7090); 29.0464 S, 115.3467 E, 0.244 km |
| Normal points | ILRS CRD v2, records 11; UTC seconds of day, two-way time of flight; epoch event 2 (ground transmit), converted to reception time |
| Files | `data/starlette_20250609.np2`, `..._20250610.np2`, `..._20250611.np2` |
| CPF | `data/starlette_cpf_250609_16101.sgf`; 180 s grid, ITRF positions |
| Span | 2025-06-09 00:17 UTC to 2025-06-11 23:33 UTC (2.74 days) |
| Source | EUROLAS Data Center (DGFI-TUM) for the International Laser Ranging Service |
| Retrieved | 2025-10-24 |

The files are the original downloads, unmodified except for the `.np2`
extension used by the EDC. Checksums (SHA-256):

    9fe2fd08df4fa91d84d5f2361adc7b440b3f6949ab6d9e67d09f9b4fc1117a73  starlette_20250609.np2
    48a384525f6a135f55d7ce86e42e5345ea76ddfe1e00d00875e8eb7e3d70b63e  starlette_20250610.np2
    32f897231cc2345ab5efd037f8e926c5a1610e45ff82771217743d6a32e06b34  starlette_20250611.np2
    cf76575d9b02490685494147fcca0c68309dee94b6e69647ed379089872616eb  starlette_cpf_250609_16101.sgf

If you use these data, acknowledge the ILRS and the Yarragadee station, and
cite:

> Pearlman, M., Noll, C., McGarry, J., et al. (2019). The ILRS: approaching
> 20 years and planning for the future. *Journal of Geodesy* 93, 2161-2180.
> https://doi.org/10.1007/s00190-019-01241-1

The CRD files also contain sessions from other stations (for example YEBL);
the parser filters on the `H2` station record and keeps YARL only.

## Method

- **CPF interpolation.** The CPF grid is 180 s. Linear interpolation is not
  adequate: over 180 s the along-track displacement is about 130 km and the
  chord/arc error is kilometres. Positions and velocities are interpolated with
  an 8-point barycentric Lagrange polynomial. The script prints the difference
  between the 8- and 10-point interpolants at grid midpoints as an
  order-of-accuracy check: the maximum is about 35 mm at interior midpoints,
  where both stencils are two-sided. At the first and last midpoints the
  8-point stencil is one-sided and the difference rises to 0.44 m; those spans
  lie outside the SLR data, which start 56 min after the first CPF node.
- **Earth orientation.** GMST from the IAU 1982 expression; UTC is used as UT1
  and polar motion is neglected. The ITRF-to-ECI velocity includes the
  Earth-rotation term `omega x r`.
- **Light time.** A laser range is half the round-trip time of flight and
  refers to the satellite at the bounce epoch. CRD record 11 in these files is
  transmit-tagged (epoch-event field 2, ground transmit), so
  `parse_slr_normal_points` converts every epoch to reception time,
  `t_rx = t_tx + tof`. The measurement model then solves the down-leg and
  up-leg light times by fixed-point iteration, keeping the satellite at the
  bounce state and rotating the station through the light time. At Starlette's
  roughly 1300 km range the bounce epoch is 4.3 ms before reception and the
  two-way range differs from the geometric range at reception by 21 m RMS.
  Applying this correction to a transmit-tagged epoch instead of converting the
  epoch first inflates the CPF residual from 7.8 m to 39.5 m; that comparison is
  how the convention was verified.
- **Force model.** J2-J5 zonal harmonics, exponential drag with Cd = 2.1 and
  A/m from the 0.12 m sphere, and SRP with Cr = 1.1, an analytic solar
  ephemeris, and a conical shadow. No tesseral harmonics, no solid Earth tides,
  no atmospheric delay, no relativity.
- **Fit.** Range-only Levenberg-Marquardt with the six-component state at the
  first normal point. `sigma_range = 10 m` sets only the weight scale; the
  state estimate is independent of it and the covariance scales as
  `sigma_range^2`. Starlette normal points are precise to centimetres, but the
  fit residuals are tens of metres because of force-model and troposphere
  systematics, so the normalized residual RMS (the square root of the reduced
  chi-square) is 2.5-5.9. The script reports both the formal covariance at the
  10 m weight and the a-posteriori covariance scaled by that factor.

## Reproduce

    julia --project=. validation/slr/run_starlette.jl

The script prints every number in `RESULTS.md` and takes a few minutes at
`reltol = abstol = 1e-11`.

## Limitations

- One satellite, one station, 2.7 days. This is a validation case, not an
  operational POD pipeline.
- The CPF is a prediction product distributed by ILRS. Its own error, including
  prediction decay after the fit span, is part of check 1 and is not separated
  from our model error.
- The station coordinates are the ILRS reference values at the epoch; station
  motion (plate tectonics) and tidal displacements are not modelled.
- The CRD epoch-event field identifies the time tag. All records used here are
  event 2 (ground transmit), so they are converted to reception time before
  comparison; the up-leg/down-leg split is modelled from the geometry rather
  than read from the file.
- Tropospheric refraction is not modelled and is not flagged as applied in
  these H4 headers. It adds about 2.4 m at zenith and 7-9 m at 15 deg
  elevation, consistent with the +4 m mean residual of check 1.
