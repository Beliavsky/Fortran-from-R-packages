# lmomco-fortran v0.1.0

A modern Fortran/FPM computational-core port of selected high-value numerical functionality from the R package **lmomco 2.5.7** by William H. Asquith.

This archive is the corrected source release. It contains actual Fortran sources in `src/`, tests in `test/`, an example in `example/`, and the original supplied R package archive in `upstream/`.

## Implemented

### Generic distribution API

```fortran
use lmomco

type(lmomco_params) :: par
real(dp) :: f, q, d

par = make_params('gev', [0.0_dp, 1.0_dp, 0.1_dp])
d = lmomco_pdf(1.2_dp, par)
f = lmomco_cdf(1.2_dp, par)
q = lmomco_quantile(0.99_dp, par)
```

The generic engine currently supports these lmomco family identifiers:

- `aep4`
- `cau`
- `emu` (finite-parameter numerical CDF/quantile)
- `exp`
- `texp`
- `gam` (two-parameter gamma)
- `gep`
- `gev`
- `gld`
- `glo`
- `gno`
- `gov`
- `gpa`
- `gum`
- `kap`
- `kmu` (finite kappa; numerical CDF/quantile)
- `kur`
- `lap`
- `lmrq`
- `ln3`
- `nor`
- `pdq3`
- `pdq4`
- `pe3`
- `ray`
- `revgum`
- `rice` (numerical CDF/quantile)
- `sla`
- `smd`
- `st3`
- `tri`
- `wak`
- `wei`

The generic API also provides RNG generation by inverse transform.

### L-moment functionality

- unbiased sample probability-weighted moments
- ordinary sample L-moments
- L-CV/L-skew/L-kurtosis conversion
- theoretical L-moments by numerical integration of a distribution quantile function
- direct L-moment fits for Normal, Exponential, Gumbel, and GEV

### Named compatibility entry points

Named PDF/CDF/quantile wrappers are supplied for the principal families:

- Exponential
- Normal
- GEV
- Gumbel
- Generalized Pareto
- Laplace
- Cauchy
- Weibull

Named L-moment routines are supplied for Exponential, Normal, GEV, Gumbel, GPA, and Weibull, together with `parexp`, `parnor`, `pargum`, and `pargev`.

### Other numerical helpers

- plotting positions
- return-period/nonexceedance transforms
- harmonic mean
- Gini mean difference
- regularized incomplete gamma and beta functions
- normal/gamma/Student-t quantiles
- fractional modified-Bessel-I series/asymptotic evaluator for fading distributions

## Scope

`lmomco` is a very large R package (hundreds of exported functions). This **v0.1.0** is a real, compilable computational-core port, but it is not a claim of one-for-one parity with every R export. In particular, specialized censored/BFR/regional estimators, many plotting/diagnostic helpers, bootstrap covariance machinery, and some specialized parameter inversions remain future targets.

The R formula/list/S3 interfaces are intentionally replaced by typed Fortran parameters and array APIs.

## Build

With FPM:

```text
fpm test
fpm run --example basic
```

The source was also validated directly using GNU Fortran:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

## Validation

The included tests cover:

- direct distribution reference identities
- GEV, GLD, and triangular CDF/quantile inversion
- sample L-moments
- theoretical Exponential L-moments
- Normal and Exponential direct L-moment fits
- named compatibility wrappers
- finite Rice and Kappa-Mu numerical evaluation

## License

The upstream `DESCRIPTION` declares `License: GPL` without a version qualifier. This port preserves that statement and includes the exact supplied package archive under `upstream/` for attribution/auditability. See `LICENSES.md`.
