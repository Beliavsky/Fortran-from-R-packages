# rfortran-core

`rfortran-core` provides small, tested modern Fortran modules for operations
that recur across translations of R packages. It is intentionally independent
of BLAS, LAPACK, operating-system APIs, and package-specific status types.

The initial modules are:

- `r_kinds`: common real and integer kinds and mathematical constants.
- `r_status`: status values returned by shared numerical procedures.
- `r_missing`: elemental IEEE missing/non-finite predicates.
- `r_descriptive`: missing-aware mean, variance, standard deviation and counts.
- `r_time_series`: biased autocovariance and autocorrelation sequences.
- `r_vectors`: R-style lagged and repeated differencing for vectors and matrices.
- `r_mod`: a compatibility facade that re-exports the initial API.

New translations should import narrowly from the defining module:

```fortran
use r_kinds, only : dp
use r_descriptive, only : r_mean, r_variance
```

The `r_mod` facade is intended to ease migration of existing translations.

## Numerical policies

- `r_mean`, `r_variance`, and `r_sd` propagate NaNs by default.
- Passing `na_rm=.true.` omits NaNs, corresponding to R's `na.rm = TRUE`.
- Empty reductions and variance with insufficient degrees of freedom return
  an IEEE quiet NaN.
- `r_variance` uses sample variance (`ddof=1`) by default.
- Autocovariances use the R `acf(..., type="covariance")` divisor `n` at every
  lag. Autocorrelations divide these values by the lag-zero autocovariance.
- Time-series routines reject NaNs and infinities and report a status rather
  than terminating the program.
- Matrix time-series operations treat rows as observations and columns as
  series. Differencing is performed along the first dimension.
