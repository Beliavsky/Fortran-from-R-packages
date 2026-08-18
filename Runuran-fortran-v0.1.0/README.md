# Runuran-fortran

Modern Fortran 2018 / FPM translation of the computational core of the R
package **Runuran 0.41**, which interfaces the **UNU.RAN** universal
non-uniform random-variate generation library.

Upstream authors: **Josef Leydold** and **Wolfgang Hoermann**.  The original
Runuran and UNU.RAN code is GPL-2.0-or-later; this translation keeps that
license and attribution.

## What is included

The translation provides native Fortran derived types for continuous,
discrete, and multivariate distributions; a reproducible xoshiro256** uniform
RNG; probability functions, quantiles and sampling; truncated distributions;
custom callback-defined distributions; mixtures; and generator objects.

Named continuous distributions include normal, beta, Cauchy, chi,
chi-square, exponential, F, Frechet, gamma, Gumbel, inverse Gaussian,
Laplace, lognormal, logistic, Lomax, Pareto, power-exponential, Rayleigh,
slash, Student t, triangular, Weibull, Burr, GIG (both parameterizations),
hyperbolic, generalized hyperbolic, variance-gamma, Meixner, and Planck.

Named discrete distributions include binomial, geometric, hypergeometric,
logarithmic-series, negative binomial, Poisson, Zipf, and arbitrary finite
probability vectors.

The user-facing UNU.RAN generator families are represented by Fortran
constructors corresponding to PINV, ARS, AROU, SROU, TDR, ITDR, TABL,
DARI, DAU, DGT, mixtures, HITRO, and VNROU.  See `TRANSLATION_NOTES.md` for
which kernels are direct algorithmic ports and which use a target-equivalent
native Fortran fallback.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic_sampling
```

No external numerical library is required.

## Minimal example

```fortran
program demo
  use runuran
  implicit none
  type(rng_state) :: rng
  type(continuous_distribution) :: d
  type(unuran_generator) :: gen
  real(dp) :: x(1000)

  call rng_seed(rng, 12345_i8)
  d = udgamma(2.5_dp, 1.2_dp)
  gen = tdr_new(d)
  call gen%sample_n(rng, x)
  print *, sum(x) / real(size(x), dp)
end program demo
```

`ud`, `up`, `uq`, and `ur` provide short Runuran-like density/PMF, CDF,
quantile, and generator-sampling entry points.  The derived types also expose
`pdf`, `logpdf`, `dpdf`, `dlogpdf`, `cdf`, `quantile`, `sample`, and
`sample_n` methods directly.
