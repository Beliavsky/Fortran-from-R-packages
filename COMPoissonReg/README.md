# COMPoissonReg-fortran

Modern Fortran 2018 translation of the computational code in the R package
`COMPoissonReg` 0.8.2, packaged for FPM.

The port covers the Conway-Maxwell-Poisson (CMP), zero-inflated CMP (ZICMP),
and zero-inflated Poisson (ZIP) probability models, the hybrid CMP
normalizing-constant algorithm, raw matrix-based regression fitting,
diagnostics, parametric bootstrap, and the package's experimental Fisher
information routines.

R formula parsing, S3 printing/summary infrastructure, Rcpp registration, data
frame handling, and vignette/plotting material are not part of the Fortran
runtime library. The complete supplied upstream package is retained under
`upstream/` for provenance.

## Build

```text
fpm build
fpm test
fpm run --example basic_cmp
```

The public umbrella module is:

```fortran
use compoissonreg
```

## Distribution example

```fortran
real(dp) :: p, mu, vr
integer :: q

p  = pcmp(4, 3.0_dp, 0.8_dp)
q  = qcmp(0.75_dp, 3.0_dp, 0.8_dp)
mu = ecmp(3.0_dp, 0.8_dp)
vr = vcmp(3.0_dp, 0.8_dp)
```

Random draws use caller-allocated integer arrays:

```fortran
integer :: x(1000)
call rcmp(size(x), 3.0_dp, 0.8_dp, x)
```

## Raw CMP regression

The R formula interface is intentionally replaced by a direct matrix API.
`lambda` and `nu` retain the package's log links.

```fortran
integer :: y(10)
real(dp) :: xmat(10,1), smat(10,1)
type(cmp_fit_t) :: fit

y = [0,1,2,1,3,2,4,1,0,2]
xmat = 1.0_dp
smat = 1.0_dp
call fit_cmp_raw(y, xmat, smat, fit)
```

For ZICMP regression, `fit_zicmp_raw` additionally accepts `W`; the
zero-inflation parameter uses the logit link. Initial values, offsets, fixed
coefficient masks, and numerical controls are represented by `cmp_init_t`,
`cmp_offset_t`, `cmp_fixed_t`, and `cmp_control_t`.

## Numerical method

The upstream package chooses between a truncated series and an asymptotic
approximation for the CMP normalizing constant. The same hybrid criterion and
default tolerances are used here. The truncation summation is performed in log
space with a recurrence between neighboring terms.

Regression uses a standalone BFGS maximizer. Its score is evaluated from the
CMP likelihood identities for `E(Y)` and `E(log(Y!))`; the covariance matrix is
based on an independent numerical Hessian at the solution. This replaces R's
`stats::optim` and `numDeriv` dependencies.

## Tests

The test suite checks:

- CMP/ZIP/ZICMP identities and limiting cases;
- independent high-precision CMP probability/moment references;
- normalizer/truncation behavior;
- random-generation moments and zero frequencies;
- raw CMP and ZICMP regression;
- residuals, AIC/BIC, equidispersion testing, leverage, deviance, and bootstrap;
- exact/Monte-Carlo-style Fisher information code paths.

See `PORTING_NOTES.md` and `API_MAP.md` for differences from the R package.
