# compositions-fortran

Modern free-format Fortran/FPM translation of the computational core of the R
package `compositions` 2.0-9.

The project focuses on numerical compositional-data analysis.  R S3 classes,
formula/model-frame dispatch, plotting, printing, and other presentation code are
not compiled.

## Build

```text
fpm build
fpm test
fpm run --example demo_compositions
```

BLAS and LAPACK are used for eigensystems, least squares, SVD/pseudoinverses,
and linear solves.  The bundled tensorA layer uses Fortran 2018 features. GNU Fortran users can
also run `scripts/validate.sh` or `scripts\validate.bat`.

## Main API

Use the facade module:

```fortran
use compositions
```

### Geometry and transforms

- `closure`, `perturb`, `power_comp`
- `clr`, `ilr`, `alr` and inverse transformations
- `apt`, `cpt`, `ipt`, `ilt` and inverses
- Helmert/basic ILR bases and balance bases
- CLR/ILR covariance and variation-matrix conversions
- Aitchison center, geometric means, pairwise log-ratios
- missing-part projection matrices

### Probability models

- Dirichlet density and simulation
- logistic-normal density and simulation
- Aitchison-distribution numerical normalizing integrals, density, and rejection simulation
- Poisson and multinomial count-composition simulation

### Estimation and multivariate statistics

- `fit_dirichlet`
- classical Aitchison covariance and PCA
- robust covariance/PCA through the supplied `robustbase` DetMCD translation
- one- and two-sample compositional normal-location likelihood-ratio statistics
- native multivariate compositional linear regression in ILR coordinates
- Mahalanobis distances and simulated empirical/max-Mahalanobis outlier calibration
- single-component outlier explanation/classification
- `PBhclust`, `PBangprox`, and source-style `PBmaxvar` principal-balance construction

### Zero and missing-data helpers

- detection-limit extraction from negative-coded values
- multiplicative-style zero replacement kernel
- MCAR simulation helper
- observed-part and summed observed-part projectors
- missing-pattern classification and conditional ALR moments
- MAR/BDL conditional imputation, projection fitting, and iterative EM-style imputation
- source-pattern imputation caches, normalized-residual Monte Carlo expectation, and source/corrected covariance outputs

### Geostatistics

- empirical log-ratio variograms
- variogram-to-log-ratio conversion
- spherical, exponential, Gaussian, linear, power, nugget, and cardinal-sine models
- complete-data ordinary compositional kriging kernel
- generalized universal compositional kriging with arbitrary trends and partially observed rows
- corrected CLR kriging covariance and source-centering compatibility mode

### Goodness-of-fit

- source-compatible Gaussian-kernel two-sample similarity statistic and permutations
- source-compatible Poisson KS statistic
- sorted-uniform simulation and Poisson-KS Monte Carlo kernel
- native energy-distance k-sample compositional test with permutation inference
- multivariate-normal energy test with parametric normal bootstrap

## Dependency translations

The full supplied translations are retained under `vendor/`.

- `robustbase-modern-fortran-0.3.0`: DetMCD is compiled and used directly.
- `bayesm-fortran-v0.1.0`: RNG/special numerical support is compiled and used directly.
- `tensorA-fortran-v0.1.0`: the complete supplied implementation is compiled into
  the library. `compositions_tensor` exposes named-axis tensors, broadcasting,
  contraction, Einstein/Riemann pairing, index dragging, batched linear algebra,
  SVD/Cholesky/pseudoinverse, and high-rank statistics.

The v0.2.0 parity pass translated the experimental missing-data and
principal-balance areas. The v0.3.0 pass adds the source-oriented imputation cache
and Monte Carlo expectation semantics, outlier calibration, energy-distance GOF,
and the full supplied tensorA numerical layer.  Because upstream
`gsiCFitWithEM` contains an empty iteration loop and `compOKriging(err=TRUE)` is
explicitly marked buggy, corrected numerical modes are provided alongside
source-compatible low-level behavior rather than reproducing undefined/broken
operations literally.

See `docs/TRANSLATION_STATUS.md` for exact coverage and the remaining source-defect/interface differences.

## Licensing

`compositions` is GPL-2-or-later.  Original source and license material are kept
under `upstream/`; dependency translations and their notices are kept under
`vendor/`.  See `LICENSES.md`.
