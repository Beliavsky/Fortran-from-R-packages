# leaps-fortran

Modern free-format Fortran translation of the computational code in the R
package `leaps` 3.2, packaged for the Fortran Package Manager (FPM).

The package performs best-subset selection for linear regression using Alan
Miller's QR updating and subset-search algorithms.  The supplied free-format
Fortran 95 versions of Miller's `lsq` and subset-search code are used as the
numerical foundation, with a new high-level Fortran API implementing the
computational behavior of the R wrapper.

## Implemented computational features

- weighted least-squares QR setup;
- intercept/no-intercept fits;
- forced-in and forced-out predictors;
- detection and handling of linearly dependent predictors;
- exhaustive leaps-and-bounds search;
- backward elimination;
- forward selection;
- sequential replacement;
- multiple best subsets per model size (`nbest`);
- RSS, R-squared, adjusted R-squared, Mallows Cp, and BIC;
- regression coefficients for a selected subset;
- covariance matrix for selected-model coefficients.

R-specific formula processing, S3 methods, printing, `biglm` object handling,
and plotting are not part of the Fortran library.  Plotting code was
explicitly omitted.

## Build and test with FPM

```text
fpm build
fpm test
fpm run --example demo_leaps
```

The source is free-form and uses explicit interfaces through modules.  The
public floating-point kind is `dp = kind(1.0d0)`.

## Basic API

```fortran
use leaps, only: dp, regsubsets_result, regsubsets_fit, get_model

type(regsubsets_result) :: fit
integer, allocatable :: ids(:)
integer :: ier
real(dp) :: rss

call regsubsets_fit(x, y, fit, nvmax=6, nbest=3, &
                    method='exhaustive', ier=ier)
call get_model(fit, 3, 1, ids, rss, ier)
```

`x` has shape `(nobs,nvar)`.  Optional `weights`, `force_in`, and `force_out`
are one-dimensional arrays of lengths `nobs`, `nvar`, and `nvar`,
respectively.  `force_in` and `force_out` are logical masks.

Methods are `exhaustive`, `backward`, `forward`, and `seqrep`.

### Model identifiers

Predictors are numbered `1:nvar`.  Identifier `0` denotes the intercept.
Thus a returned model `[0,1,3]` means intercept plus predictors 1 and 3.
The `size_predictors` argument to `get_model` and `model_coefficients`
**excludes** the intercept.

### Result arrays

The first dimension of `fit%rss`, `fit%rsq`, `fit%adjr2`, `fit%cp`,
`fit%bic`, and `fit%valid` is the total number of fitted columns, including
the intercept when present.  The second dimension is the rank from 1 through
`nbest`.  Use `fit%valid` before reading a slot.

`fit%model(1:p,p,rank)` stores the model identifiers for a valid model with
`p` fitted columns.

## Error status

`regsubsets_fit` reports zero on success.  Positive values denote invalid
arguments or setup/search failures.  The result object also stores the last
status in `fit%status`.

## Numerical provenance

- `src/leaps_lsq.f90`: modernized package-local copy of Alan Miller's
  free-format least-squares module.
- `src/leaps_find_subsets.f90`: modernized package-local copy of Alan Miller's
  free-format subset-search module.
- `src/leaps.f90`: high-level API translating the computational R wrapper.
- `upstream/`: verbatim supplied Miller files plus the relevant original
  `leaps` R/fixed-form Fortran sources for comparison and attribution.

See `LICENSE.md` for licensing and provenance details.
