# Porting notes

## Source baseline

- R package: `distr`
- Upstream version: 2.9.7
- Upstream package date: 2025-01-11
- Uploaded source archive SHA-256: `9c68fbfdc1c7068dd4387ce03978afffc97e0f0d17d8c96f2e9b590d46e2b853`
- Upstream declared package license: LGPL-3
- Upstream `src/ks.c`: GPL-2.0-or-later (R Core-derived)

The original DESCRIPTION, NAMESPACE, CITATION, KS source, QQ sources, and the R sources most directly relevant to convolution, lattice distributions, and mathematical transformations are retained under `upstream/`.

## Design mapping

R `distr` represents distributions using S4 classes whose slots contain `r`, `d`, `p`, and `q` functions. Fortran has no direct equivalent to that dynamic S4 model, so this port uses a recursive value type, `distribution_t`. Each value records a named distribution or a composition tree. Type-bound methods evaluate density/mass, CDF, survival probability, quantile, random variates, and moments.

This permits native expressions such as

```fortran
z = normal_dist(0.0_dp, 1.0_dp) + exponential_dist(1.0_dp)
y = truncate_dist(z, lower=-2.0_dp, upper=3.0_dp)
```

without an R interpreter or procedure-closure runtime.

## Named distributions and special functions

Native PMF/PDF, CDF, quantile, RNG, and common moment behavior are provided for Dirac, binomial, hypergeometric, Poisson, negative binomial, geometric, uniform, normal, lognormal, Cauchy, F, Student t, chi-square, exponential, Laplace, gamma, beta, logistic, Weibull, and arcsine distributions.

Noncentral beta, F, t, and chi-square are implemented without Rmath. Noncentral beta and chi-square use Poisson mixtures; both lower and upper tails are summed directly in v0.2.0. Noncentral F is reduced to noncentral beta. Noncentral t uses deterministic quadrature over the defining noncentral-normal/chi-square representation, with its survival probability evaluated through the reflected noncentral-t CDF identity.

The port includes the incomplete beta/gamma functions required for CDFs; inverse normal; digamma, trigamma and inverse digamma; and matrix square-root/generalized-inverse/linear-solve helpers corresponding to the computational portions of the upstream package.

## Distribution construction and algebra

Translated concepts include `DiscreteDistribution`, `EmpiricalDistribution`, weighted empirical laws, computational lattice distributions, mixtures, mixed discrete/continuous distributions, affine transforms, `Truncate`, `Minimum`, `Maximum`, convolution, `convpow`, `CompoundDistribution`, `Huberize`, `abs`, `exp`, `log`, square root, reciprocal, and powers.

`lattice_dist(offset,step,prob)` is a compact constructor over the same exact discrete-support representation; the full R S4 lattice metadata hierarchy is intentionally not duplicated.

Compatible operations are simplified before a generic composition node is built. Version 0.2.0 recognizes additional closed forms for Dirac shifts, normal, Poisson, Cauchy, binomial with common probability, negative binomial with common probability, gamma with common scale, equal-rate exponential, compatible gamma-plus-exponential, and central chi-square convolution. Affine transforms of several named families are similarly collapsed to named distributions when mathematically exact.

## FFT convolution in v0.2.0

Upstream `distr` contains a general-purpose FFT convolution implementation. Version 0.2.0 adds a native dependency-free radix-2 complex FFT in `distr_fft.f90` and uses the convolution theorem in two paths:

1. `convolve_fft(a,b,grid_points,tail_prob)` discretizes continuous probability mass over finite quantile bounds on a common cell width, convolves the masses by FFT, and returns a normalized piecewise-uniform grid distribution.
2. Compatible discrete lattice supports are convolved directly by FFT on their PMF arrays. Infinite lattice laws are truncated according to `tail_prob` before the FFT.

`convpow_fft` raises one Fourier transform to the requested convolution power. `convpow(base,n,grid_points=...)` now uses the FFT path for non-analytic continuous powers; in v0.1.0 those optional arguments were retained only for API correspondence and ignored.

The FFT discretization is not intended to reproduce every internal interpolation detail of upstream `distr::convpow`. It implements the same mathematical FFT-convolution strategy in native Fortran. `convolve` remains available as the accuracy-oriented adaptive-integration path for generic continuous sums.

## Tail and log-probability behavior in v0.2.0

Version 0.1.0 supported R-like `lower_tail`/`log_p` switches but obtained an upper tail as `1-cdf`. Version 0.2.0 adds direct methods

```fortran
p  = d%sf(x)
lp = d%logcdf(x)
ls = d%logsf(x)
```

and direct upper-tail formulas for the named laws. This avoids cancellation such as `1.0_dp-normal_cdf(10)`, which is zero in binary64 even though the mathematical survival probability is nonzero.

Upper-tail quantiles are also solved directly. For example,

```fortran
q = d%quantile(log(1.0e-50_dp), lower_tail=.false., log_p=.true.)
```

does not form `1-1e-50`.

Major named laws also have direct log-density formulas, so a request for a log density need not first underflow the ordinary density to zero.

## Generic transformations in v0.2.0

The original exact `exp` and `abs` transforms are supplemented by:

- `log_transform`;
- `sqrt_transform`;
- `reciprocal_transform`;
- `power_transform`.

Positive-domain powers permit arbitrary real exponents. Positive integer odd powers use one signed inverse branch; positive even powers use both inverse branches. Negative powers require positive support with no atom at zero. These restrictions deliberately avoid silently constructing complex-valued or singular transformations in a real univariate distribution type.

The implementation propagates PMF/PDF, CDF, survival probability, quantiles, RNG behavior, and support bounds through the supported transformations.

## Mixed distributions and left limits

The upstream package distinguishes Lebesgue components and atoms explicitly. The Fortran value representation is lighter, but v0.2.0 improves `cdf_left` for mixtures, affine transforms, truncation, min/max, transformations, Huberization, and convolutions involving discrete components. This is important for upper-tail and transformation formulas at atoms.

## Generic moments

In addition to mean, variance, and standard deviation, v0.2.0 adds raw moments, central moments, skewness, and excess kurtosis. First and second moments use existing analytic formulas where available. Higher moments use exact finite/truncated support sums for discrete laws and adaptive integration over small omitted tails for continuous laws.

## Kolmogorov-Smirnov and QQ confidence bands

The numerical algorithms from upstream `src/ks.c` are translated to `distr_ks.f90`: asymptotic two-sided KS probability, exact two-sample Smirnov probability, and exact one-sample Kolmogorov probability. This module retains the upstream GPL-2.0-or-later license.

The non-graphical computations from exported `qqbounds` and its helpers are translated to `distr_qq.f90`: exact/asymptotic KS critical values, exact symmetric and shortest-asymmetric binomial pointwise intervals, asymptotic pointwise intervals, and simultaneous/pointwise QQ bounds. Upstream's asymptotic `.q2kolmogorov` complement convention is intentionally reproduced for behavioral compatibility. The QQ module is GPL-3.0-or-later because it depends on `distr_ks`.

## Deliberately omitted R-only areas

The following are not native numerical algorithms and remain omitted:

- S4 class/slot registration, coercion, and replacement-method plumbing;
- printing, `show`, display labels beyond a simple Fortran `label`, and startup messages;
- `plot`, QQ-plot rendering, graphical device management, and graphics helpers;
- package option management and old-object version conversion;
- namespace masking helpers and R environment/introspection utilities.

## Remaining computational differences

Important gaps after v0.2.0 are:

1. **Arbitrary persistent user closures.** R can store arbitrary user `d/p/q/r` functions in distribution objects. The Fortran value type currently supplies named analytic laws plus discrete/grid/KDE constructors rather than persistent procedure closures.
2. **The complete R `Math` generic.** Log, exp, abs, sqrt, reciprocal, and useful real powers are available, but many-to-one periodic functions and special-function transformations such as `sin`, `gamma`, `lgamma`, and `digamma` are not represented as general distribution transforms.
3. **Exact upstream FFT interpolation details.** The v0.2 FFT engine follows the same convolution-theorem strategy but has its own cell-mass discretization and interpolation representation.
4. **Full S4 representation utilities.** `decomposePM`, `flat`, `simplifyD`, old-version conversion, symmetry/space metadata, and similar object-management facilities are not duplicated.
5. **Parallel execution.** No OpenMP layer or package-wide thread setting is introduced; the library remains serial and dependency-free.

## Validation

Four test programs cover named central/noncentral distributions, composition and arithmetic, RNG moments, KS probabilities, QQ bounds, mixed-distribution left limits, extreme survival probabilities, direct upper-tail log quantiles, log densities, transformations, closed-form convolution simplification, continuous FFT convolution, lattice FFT convolution, FFT convolution powers, and generic moments.

The release was compiled and all tests run successfully with GNU Fortran 14.2 using Fortran 2018, warnings, implicit-interface diagnostics, bounds/runtime checking, and backtraces. The final strict validation does not require a nonstandard free-line-length option; all Fortran source lines are at most 132 columns.

FPM itself was not installed in the translation environment, so the FPM manifest/layout was validated and the same `src/`, `test/`, and `example/` sources were compiled and linked directly with gfortran.
