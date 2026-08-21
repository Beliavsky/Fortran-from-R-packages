# Validation

The package is validated directly with GNU Fortran using:

```text
-std=f2018
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-Wimplicit-interface
-Werror=implicit-interface
```

Fortran 2018 is used because the compiled supplied tensorA implementation uses
features beyond the earlier Fortran-2008 validation profile.

The combined suite contains 19 independent test programs plus the demo. It
checks:

- closure and transform round trips
- CLR/ILR/variation covariance identities
- PBhclust, PBangprox, and source PBmaxvar basis orthogonality
- Dirichlet density/RNG/fitting and logistic-normal calculations
- classical and DetMCD robust compositional statistics
- exact ILR regression recovery
- zero/detection-limit and missing projector helpers
- MAR/BDL conditional ALR imputation and projection/iterative fitting
- source-pattern imputation caches, normalized-residual Monte Carlo acceptance,
  and source versus corrected covariance conventions
- empirical/max Mahalanobis calibration and outlier classification
- variograms and complete/missing-component kriging interpolation/covariance
- registered kernel-density and Poisson-KS routines
- energy-distance k-sample and multivariate-normal GOF tests
- rank-3 named tensor bridging
- all six supplied tensorA regression programs covering named algebra, geometry,
  linear algebra, batched statistics, and miscellaneous high-rank operations

FPM was not available in the translation environment, so the FPM project layout
was validated by compiling the same `src`, `test`, and `example` files directly
with GNU Fortran and linking BLAS/LAPACK.
