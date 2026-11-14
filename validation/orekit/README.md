# Independent propagation cross-validation with Orekit

The reference states in `reference/` were generated with **Orekit 13.1.9**
(Apache-2.0) running on OpenJDK 21, not with this package. `OrekitReference.java`
propagates the same initial state and force model as the benchmark scenario:

| Item | Value |
|---|---|
| Initial state | `KeplerianState(6900 km, 0.002, 52 deg, 30 deg, 40 deg, 0)` in EME2000 |
| Epoch | 2000-01-01T12:00:00 TAI |
| Force models | Newtonian + `J2OnlyPerturbation`, optionally `DragForce` |
| J2 | 1.0826269e-3, equatorial radius 6378.137 km, mu 398600.4418 km^3/s^2 |
| Drag | spherical Earth (no flattening), `SimpleExponentialAtmosphere(2.789e-10 kg/m^3 at 200 km, H = 50 km)`, Cd 2.2, area/mass 0.01 m^2/kg |
| Integrator | Dormand-Prince 8(5,3), absTol 1e-6 m, relTol 1e-13 |
| Output grid | 0 to 86400 s, every 600 s |

The Julia side uses `Tsit5` with `reltol = abstol = 1e-12`. `compare.jl`
reproduces the comparison from the repository root:

    julia --project=. validation/orekit/compare.jl

Measured maximum differences over one day:

| Case | max position | max velocity |
|---|---|---|
| J2 | 0.0012 m | 0.0013 mm/s |
| J2 + drag | 0.0028 m | 0.0031 mm/s |

The test suite asserts a looser bound (5 cm, 0.01 mm/s) so that solver and
platform variations do not make CI brittle.

## Reproducing the reference

    # jars (Maven Central)
    mkdir -p lib && cd lib
    for a in orekit-13.1.9 hipparchus-core-4.0.1 hipparchus-geometry-4.0.1 \
             hipparchus-ode-4.0.1 hipparchus-fitting-4.0.1 hipparchus-optim-4.0.1; do
      case $a in orekit*) g=org/orekit/orekit/13.1.9 ;;
                   hipparchus*) g=org/hipparchus/${a%-4.0.1}/4.0.1 ;; esac
      curl -LO "https://repo1.maven.org/maven2/$g/$a.jar"
    done
    cd ..

    # Orekit data (Earth orientation, leap seconds); not needed at test time
    curl -L -o orekit-data-main.zip \
      https://gitlab.orekit.org/orekit/orekit-data/-/archive/main/orekit-data-main.zip
    unzip -q orekit-data-main.zip

    javac -cp "lib/*" OrekitReference.java
    java -cp "lib/*;." OrekitReference orekit-data-main j2     reference/reference_j2.csv
    java -cp "lib/*;." OrekitReference orekit-data-main j2drag reference/reference_j2drag.csv

Reference file checksums (SHA-256):

    bc3e0ac34fa3ea092573f3149ac921932863b7c1438814166849510656eb8920  reference/reference_j2.csv
    d6ff07917642370ad392158684d4cedf9a82eea3a8d416b159f9ec1800693283  reference/reference_j2drag.csv

## Scope

This validates the numerical propagation and the J2 and drag force models
against an independent implementation. It does not validate the solar radiation
pressure model, the observation model, or the estimator. The committed CSVs are
the reference of record; the Orekit jars and data directory are not committed
and are only needed to regenerate them.

Orekit is distributed under the Apache License 2.0. The Orekit data directory is
a collection of publicly available IERS, JPL and IAU files; see the
`orekit-data` project for the individual terms.
