# Translation notes

## Source package

- R package: CEoptim
- Version: 1.3
- Date: 2023-10-04
- License declared by package: GPL (>= 2.0)

The package has `NeedsCompilation: no`; all computational code is written in R.
The original sources used for this translation are retained under `original/`.

## File mapping

### `R/CEoptim.R`

Translated to `src/ceoptim_core.f90`, with public data structures in
`src/ceoptim_types.f90`.

Preserved computational behavior includes:

- continuous, discrete, and mixed CE sampling;
- initial elite fit before entering the main iteration loop;
- unbiased sample standard deviation, matching R's `sd`;
- independent Gaussian covariance `diag(sd**2)`;
- smoothing formulas for mean, standard deviation, and category probability;
- best-so-far optimizer/value retention;
- `gammat` tracking as implemented by CEoptim;
- SD/probability convergence thresholds;
- the package's no-improvement stopping logic;
- 0-based categorical values;
- maximization by sign reversal;
- iteration state/probability storage.

The R package reports `nfe=iter*N`, which omits the initial sample evaluation.
`ce_result%nfe` preserves that behavior. `ce_result%actual_nfe` additionally
reports the actual callback count, `(niter+1)*N` for a normal completed run.

### `R/rtmvnorm.R`

Translated to `src/ceoptim_sampling.f90`.

The two-stage algorithm is preserved:

1. multivariate-normal accept/reject with adaptive proposal batch size;
2. Gibbs fallback when acceptance is inadequate.

Differences that do not change the target distribution:

- covariance square roots and pseudoinverses use LAPACK symmetric
  eigendecomposition rather than `MASS::mvrnorm`/`MASS::ginv`;
- the Gibbs conditional mean/variance use the equivalent precision-matrix
  formula directly;
- univariate truncated normals use inverse-CDF sampling in the central region
  and Robert-style exponential rejection in extreme tails instead of
  `msm::rtnorm`.

The public sampler accepts a covariance matrix, as the R function does.
Positive-definite covariance is the normal Gibbs use case; semidefinite
covariance can be sampled in the accept/reject stage.

### `R/dirichletrnd.R`

Translated directly: independent gamma(shape=`a(j)`, rate=1) draws are
normalized by their row sum. Gamma variates use Marsaglia-Tsang sampling.

### `R/print.CEoptim.R`

Not translated. It is presentation/S3 infrastructure rather than numerical
code. Fortran callers access fields of `ce_result` directly.

## R dependencies removed

- `MASS::mvrnorm` -> native covariance factor + normal RNG
- `MASS::ginv` -> LAPACK symmetric pseudoinverse
- `msm::rtnorm` -> native truncated-normal sampler
- `stats::rgamma` -> native gamma RNG
- `stats::sd` -> explicit unbiased SD calculation
- `sna::gplot` -> unused by the numerical core and omitted

## Random streams

R and Fortran intentionally use different RNG engines. Algorithmic regression
is therefore based on convergence, feasibility, distribution identities, and
known objective values rather than byte-for-byte sample equality.

## Validation flags

The release is tested with GNU Fortran using both an optimized build and a
bounds-checked warning-heavy build including:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

The test programs use internal objective callbacks; GNU ld may consequently
print an executable-stack warning for those test executables. This comes from
the compiler's callback trampoline and is not emitted by the library objects.
