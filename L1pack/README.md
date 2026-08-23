# L1pack-fortran

Modern Fortran translation of the numerical core of the R package **L1pack** 0.62-4.

The project uses free-form Fortran 2018 and the Fortran Package Manager (FPM). The
`fastmatrix-fortran` v0.2.0 translation is vendored as a local dependency because
upstream L1pack directly depends on fastmatrix.

## Implemented numerical API

- Univariate Laplace: `dlaplace`, `plaplace`, `qlaplace`, `rlaplace`
- Multivariate Laplace: `dmlaplace`, `log_dmlaplace`, `rmlaplace`
- Algorithm 478 / Barrodale-Roberts L1 regression: `l1fit`, `lad_fit_br`
- Phillips IRLS/EM LAD regression: `lad_fit_em`
- LAD dispatcher and utilities: `lad_fit`, `vcov_lad`, `confint_lad`,
  `predict_lad`, `simulate_lad`, `lad_quantile_residuals`, `lad_deviance`
- Multivariate Laplace EM fitting: `laplace_fit`
- Equal-means restricted multivariate Laplace fit: `laplace_fit_equal`
- Generalized spatial median / Kotz-type fit: `spatial_median_fit`
- L1 concordance correlation: `l1ccc`, `laplace_rho1`, `gaussian_rho1`,
  `ustat_rho1`, `l1ccc_bootstrap`
- Laplace Wilson-Hilferty transform: `wh_laplace`
- Non-plotting simulated diagnostic envelope: `envelope_laplace`

The main aggregate module is:

```fortran
use l1pack
```

## Build

```text
fpm test
fpm run --example demo
```

LAPACK and BLAS are linked because the vendored fastmatrix compatibility layer
uses them.

A direct strict GNU Fortran build was also used during translation with:

```text
-std=f2018 -Wall -Wextra -Werror -fcheck=all
```

## Example

```fortran
use l1pack

real(dp) :: x(5,2), y(5)
type(lad_result) :: fit

x(:,1) = 1.0_dp
x(:,2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
y = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 100.0_dp]

call lad_fit_br(x, y, fit)
print *, fit%coefficients
```

The LAD solution is approximately `(1, 2)`, despite the final outlier.

## Scope

R formula/model-frame handling, S3 methods, console formatting, and plotting are
not reproduced. Their numerical work is exposed directly through Fortran arrays
and derived result types. `envelope.Laplace` is retained as a numerical envelope
routine but does not draw a plot.

See `API_MAP.md` and `PORTING_NOTES.md` for details.
