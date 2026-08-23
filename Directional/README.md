# Directional-fortran

Modern Fortran/FPM port of the numerical core of CRAN **Directional 7.7** (2026-07-20).
The upstream package declares `GPL (>= 2)`, so this port is distributed as GPL-2.0-or-later.

## Implemented through v0.4.0

- Geometry: latitude/longitude <-> unit Cartesian coordinates, geodesic rotation matrices, Haversine distance matrices, row normalization.
- Circular distributions: von Mises, multimodal von Mises, cardioid, wrapped Cauchy, wrapped normal, circular beta, circular exponential, circular Purkayastha; wrapped-Cauchy CDF and numerical von-Mises CDF.
- Spherical distributions: von Mises-Fisher, spherical Cauchy, Poisson-kernel based distribution, Purkayastha, isotropic angular Gaussian, 3-D ESAG density plus arbitrary-dimensional ESAG density and simulation.
- Random generation: von Mises-Fisher, von Mises, spherical Cauchy, PKBD, arbitrary-dimensional ESAG, and spherical-Cauchy/PKBD mixtures.
- Statistics: circular mean/resultant summaries and kappa estimate, circular-circular correlation coefficient, distance correlation (generic, circular, spherical), spherical association diagnostic, and sample spherical median (medoid definition).
- FPM tests for normalization, uniform special cases, coordinate transforms, and rotations.

Plotting, `rgl`, `sf`, Natural Earth map wrappers, ggplot contours, and R parallel/foreach orchestration are intentionally omitted.

## Build

```text
fpm test
```

## Dependency translations supplied with the porting task

Translations of Rfast, Rfast2, rangen, Rnanoflann, and bigstatsr were inspected. This first Directional port keeps its public numerical core standalone so that it compiles without forcing all dependency packages into one build graph. `DEPENDENCY_MAP.md` records where their APIs correspond to upstream calls and which later parity targets can delegate to them.


## v0.2 numerical parity

Version 0.2 adds inference, mixture, classification, uniformity-test, Fisher-Bingham/Kent, Bingham, and matrix-Fisher modules. See `API_MAP.md` and `PORTING_NOTES.md` for exact parity notes and remaining targets.

## v0.3 numerical parity

Version 0.3 adds general d-dimensional ESAG density/simulation, PKBD simulation, spherical-Cauchy and PKBD mixture density/simulation/EM fitting, Kent MLE, generic/circular/spherical distance correlation, and multivariate KNN regression with cross-validated tuning. `test_parity_v03` exercises these additions under runtime bounds checking.

## v0.4 numerical parity

Version 0.4 adds general isotropic angular Gaussian MLE, spherical isotropic projected-Cauchy MLE, circular isotropic projected-Cauchy MLE, GCPC MLE, and reusable embedding/high-concentration/heterogeneous two-sample permutation tests. `test_parity_v04` exercises these additions under runtime bounds checking.

## v0.5.0

v0.5.0 closes the main remaining computational parity gaps: IAG/SIPC/CIPC/spherical-Cauchy/PKBD/SPML/GCPC/ESAG/SESPC regression, spherical-spherical rotation regression, general-group embedding/high-concentration/heterogeneous/likelihood-ratio ANOVA, and centered two-sample bootstrap tests. The complete accumulated test suite passes with gfortran runtime checking.

For regression, pass the design matrix explicitly; add a leading column of ones if an intercept is desired.
