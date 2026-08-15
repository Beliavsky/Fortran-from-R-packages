# distr-fortran

`distr-fortran` is a modern Fortran/FPM computational port of version 2.9.7 of the R package **distr** (Object Oriented Implementation of Distributions).

The port removes the R runtime, S4 class machinery, graphics, printing, and startup infrastructure and keeps the reusable probability/statistics algorithms in a native Fortran API. Version 0.2.0 adds FFT convolution, a broader transformation algebra, and tail-safe probability methods.

## Build

```text
fpm build
fpm test
fpm run --example basic
fpm run --example v02_features
```

The sources use free-form Fortran 2018 and have no external runtime dependencies.

## Basic API

```fortran
use distr

type(distribution_t) :: x, y, z
real(dp), allocatable :: sample(:)

x = normal_dist(mean=1.0_dp, sd=2.0_dp)
y = exponential_dist(rate=0.5_dp)
z = x + y

print *, x%density(0.0_dp)
print *, x%cdf(0.0_dp)
print *, x%sf(8.0_dp)
print *, x%logsf(8.0_dp)
print *, x%quantile(0.95_dp)
print *, z%mean(), z%variance()

call seed_rng(12345)
sample = z%random(1000)
```

`density`, `cdf`, and `quantile` retain R-like log/lower-tail options. The dedicated `sf`, `logcdf`, and `logsf` methods avoid forming `1-cdf` in small upper tails. Upper-tail quantiles are solved directly, including when the supplied probability is logarithmic.

Vector convenience routines are provided for density, CDF, survival probability, log-CDF, log-survival, and quantiles.

## Named distributions

The following R `distr` families have native constructors:

| R constructor | Fortran constructor |
|---|---|
| `Dirac` | `dirac_dist` |
| `Binom` | `binomial_dist` |
| `Hyper` | `hypergeometric_dist` |
| `Pois` | `poisson_dist` |
| `Nbinom` | `negative_binomial_dist` |
| `Geom` | `geometric_dist` |
| `Unif` | `uniform_dist` |
| `Norm` | `normal_dist` |
| `Lnorm` | `lognormal_dist` |
| `Cauchy` | `cauchy_dist` |
| `Fd` | `f_dist` |
| `Td` | `student_t_dist` |
| `Chisq` | `chisq_dist` |
| `Exp` | `exponential_dist` |
| `DExp` | `laplace_dist` |
| `Gammad` | `gamma_dist` |
| `Beta` | `beta_dist` |
| `Logis` | `logistic_dist` |
| `Weibull` | `weibull_dist` |
| `Arcsine` | `arcsine_dist` |

Noncentral beta, F, t, and chi-square parameters are supported without Rmath.

## Distribution composition

The library includes:

- arbitrary finite discrete, lattice, empirical, and weighted-empirical distributions;
- finite mixtures and mixed discrete/continuous distributions;
- affine transformations and overloaded `+`, `-`, `*`, `/` with scalars;
- truncation, minimum, maximum, convolution, convolution powers, and compound distributions;
- `exp`, `log`, `sqrt`, reciprocal, absolute-value, and general power transformations where the real-valued inverse branches are well-defined;
- Huberization/winsorization;
- KDE/grid distributions built from simulations.

Compatible closed-form sums are simplified automatically. This includes Dirac shifts, normal, Poisson, Cauchy, common-probability binomial and negative-binomial, common-scale gamma, equal-rate exponential, gamma-plus-compatible-exponential, and central chi-square sums. Compatible affine transformations are also kept in named closed form when possible.

## FFT convolution

`convolve` remains the accuracy-oriented generic operation: discrete laws use support summation and continuous laws use adaptive integration if no closed form is available.

Version 0.2.0 adds a dependency-free radix-2 FFT engine:

```fortran
z = convolve_fft(x, y, grid_points=4096, tail_prob=1.0e-10_dp)
w = convpow_fft(x, 20, grid_points=4096)
```

For continuous laws, probability masses are discretized on a common finite quantile grid, convolved with the FFT, and returned as a piecewise-uniform grid distribution. For compatible lattice laws, mass convolution uses the FFT directly. `convpow(base,n,grid_points=...)` now honors its grid argument and selects the FFT path for non-analytic continuous convolution powers.

## Transformations

```fortran
x = log_transform(lognormal_dist(0.0_dp, 1.0_dp))
y = sqrt_transform(gamma_dist(3.0_dp, 2.0_dp))
z = power_transform(normal_dist(), 2.0_dp)
w = reciprocal_transform(exponential_dist())
```

Positive-domain powers may use arbitrary real exponents. Positive integer odd/even powers also support signed input with the appropriate one- or two-branch inverse. Negative powers require positive support with no atom at zero.

## Moments

In addition to `mean`, `variance`, and `sd`, distributions expose:

```fortran
m3 = x%raw_moment(3)
c4 = x%central_moment(4)
sk = x%skewness()
ku = x%excess_kurtosis()
```

Closed-form first and second moments are used where available; higher generic moments use support sums or adaptive integration.

## Other numerical code

`distr_special` implements incomplete beta/gamma functions, inverse normal, digamma/trigamma/inverse-digamma calculations, and noncentral distribution calculations.

`distr_matrix` contains positive-semidefinite matrix square-root/generalized-inverse helpers.

`distr_ks` contains translated Kolmogorov-Smirnov probability routines originating in R Core and has a separate GPL license boundary. Import it explicitly:

```fortran
use distr_ks, only : p_ks2_asymptotic, p_smirnov2x, p_kolmogorov2x
```

The non-graphical QQ confidence-band calculations are in `distr_qq`; plotting itself is omitted.

## Scope

This is a computational translation, not an implementation of R's S4 object system. Plotting, QQ rendering, printing/show methods, startup messages, compatibility/version conversion, R option management, and similar R-only infrastructure are deliberately omitted. Arbitrary persistent R-style procedure closures and the entire R `Math` generic are not reproduced; see `PORTING_NOTES.md` for details.

## Licensing and attribution

The `distr`-derived Fortran modules are licensed under LGPL-3.0-only, following the upstream package declaration. `src/distr_ks.f90` is separately GPL-2.0-or-later because it translates `distr/src/ks.c`, itself taken from R Core; `src/distr_qq.f90` is GPL-3.0-or-later because it depends on that code. See `LICENSES.md`, the `COPYING*` files, and `upstream/` for provenance and citation information.
