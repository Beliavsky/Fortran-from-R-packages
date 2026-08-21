# Translation status

Source: R package `compositions` 2.0-9.

## Implemented numerical areas

### Compositional geometry

Implemented: closure, perturbation, powering, CLR, ILR, ALR, APT, CPT, IPT,
ILT and inverses; ILR/balance bases; covariance-coordinate conversions;
variation matrices; Aitchison/geometric means; pairwise log-ratios.

All three principal-balance choices exposed by `gsi.PrinBal` have native
Fortran algorithms: Ward/variation-distance `PBhclust`, exhaustive angular
`PBangprox`, and the source `PBmaxvar` SVD/sign-refinement state machine.

### Distributions, estimation, and statistics

Implemented: Dirichlet, logistic-normal, and Aitchison models; Poisson and
multinomial count compositions; Dirichlet Newton fitting; classical and DetMCD
robust covariance/PCA; normal-location tests; ILR multivariate regression;
prediction; Mahalanobis distances.

The outlier-calibration layer includes simulated empirical and maximum
Mahalanobis distributions, calibrated quantiles/CDFs, corrected/uncorrected
outlier decisions, and single-component explanation classification. Robust
calibration uses the bundled DetMCD implementation.

### Zero, missing, censored, and imputation computation

Implemented: detection-limit helpers, zero replacement, MCAR generation,
missing-pattern classification/indexing, observed-part projectors, conditional
ALR moments, MAR/BDL conditional imputation, projection fitting, and iterative
regression/imputation.

The v0.3 source-oriented cache path additionally implements the internal
pattern structure used by `gsiCImpAcompClrExpectation`: missing-order maps,
conditional predictor matrices, residual covariance/Cholesky factors, supplied
normalized-residual Monte Carlo draws, BDL acceptance, accepted-draw counts,
and pattern reuse. It returns the covariance accumulated in the source ALR
coordinate convention and, separately, a corrected CLR-projected covariance.

Two upstream experimental-source defects are not reproduced as undefined
memory behavior. The projection C routine computes a design entry but omits it
from the multiplication; a `source_compatible` mode remains available while
the default uses the intended design-weighted calculation. In the no-accepted
Monte Carlo branch of `gsiCImpAcompClrExpectation`, the source indexes beyond
the allocated ALR vector; the Fortran version uses the bounded conditional mean
and corresponding finite outer-product covariance instead. The upstream
`gsiCFitWithEM` iteration loop is empty in 2.0-9, so `fit_acomp_em` is a useful
completed extension rather than a claim to reproduce a nonexistent iteration.

### Geostatistics

Implemented: log-ratio variograms, variogram conversions/models, complete-data
ordinary compositional kriging, and generalized `gsiCGSkriging`-style
universal kriging with arbitrary trends and row-specific ALR references for
partially observed compositions. Only finite, strictly positive components
enter each row's ALR system. A corrected CLR kriging covariance is provided;
a compatibility option retains the source final-centering convention without
reproducing its documented invalid indexing.

### Goodness-of-fit

Implemented the package's registered native kernels: `gsiDensityCheck`,
`gsiKSPoisson`, sorted-uniform simulation, and Poisson-KS simulation, including
the original lower-triangle normalization of `gsiDensityCheck`.

Also implemented native numerical counterparts for the external energy-package
wrappers: k-sample energy distance with permutation inference and multivariate
normal energy testing after covariance standardization with a parametric normal
bootstrap. The normal-reference expectation is evaluated through an equivalent
noncentral-chi Poisson mixture, avoiding a GSL dependency.

### High-rank tensor algebra

The complete supplied `tensorA-fortran-v0.1.0` implementation is compiled into
the main library. Through `compositions_tensor`, users have dynamic named-axis
tensors, axis lookup/reorder/rename, named contraction and broadcasting,
Einstein/Riemann pairing, index dragging, inverse/solve/SVD/Cholesky/power,
pseudoinverse/norms, and high-rank mean/variance/covariance operations. The six
supplied tensorA regression programs are part of the combined test suite.

## Dependency use

The supplied `robustbase` DetMCD and `bayesm` RNG translations are compiled into
the package. The complete supplied `tensorA` translation is likewise compiled
and is no longer only a vendored provenance dependency.

## Remaining differences

No additional public numerical method from the previously identified gap list
is intentionally omitted. Remaining differences are implementation/interface
level rather than missing statistical algorithms:

1. Undefined or unfinished experimental upstream C behavior is represented by
   explicit source-compatible or corrected finite modes rather than undefined
   memory accesses or an empty top-level EM loop.
2. The energy tests use the same statistics/bootstrap structures but a native
   RNG and an algebraically equivalent normal-reference expectation instead of
   R/energy's RNG and GSL hypergeometric evaluator; tiny Monte Carlo/floating
   differences are therefore expected under different seeds/platforms.
3. R session caches, exact RNG streams, names/dimnames/factor/date metadata,
   recycling, `NA` payload distinctions, S3/formula/model-frame dispatch, and
   plotting/report infrastructure are not emulated by the numerical library.

## Excluded non-numerical infrastructure

R S3 classes (`acomp`, `rcomp`, `aplus`, `rplus`, `ccomp`, `rmult`), formula
processing, data frames, plotting/biplots, printing, interactive menus, and
exact R object metadata semantics are omitted. Their numerical transforms and
estimators are exposed through typed arrays and explicit missing-type inputs.
