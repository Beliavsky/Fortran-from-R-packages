# Translation notes

## Upstream

- R package: Runuran 0.41 (2025-04-07)
- Bundled engine: stripped UNU.RAN sources, SVN revision 8030
- Authors: Josef Leydold and Wolfgang Hoermann
- License: GPL version 2 or later

The upstream archive contains roughly 147 C files in `src/unuran-src` plus the
R/S4 interface.  This translation targets the computational behavior rather
than reproducing the R object system, `.Call` glue, parser, serialization, or
printing machinery.

## Native Fortran numerical layer

The translation is self-contained.  It includes regularized incomplete beta
and gamma functions, central beta/gamma/chi-square/t/F probability functions
and quantiles, adaptive quadrature, root finding, modified Bessel K for real
order, and a Lanczos complex-log-gamma magnitude used by the Meixner density.
No Rmath, GSL, BLAS, LAPACK, or C bridge is required.

A xoshiro256** generator with SplitMix64 seeding supplies the base uniform
stream.  Its modulo-2^64 additions and multiplications are implemented with
Fortran bit intrinsics, avoiding reliance on undefined signed-integer overflow.  Standard normal, gamma, beta, Poisson, binomial, negative-binomial,
chi-square, t, F, Cauchy and inverse-Gaussian generators are implemented
natively.

## Generator-method mapping

- **ARS**: native adaptive rejection sampling with tangent envelopes for
  log-concave univariate densities.  Numerical derivatives are used when an
  analytic callback is not supplied.
- **PINV/HINV/NINV**: native numerical CDF inversion.  This preserves the
  target distribution to the numerical integration/root-finding tolerance,
  but does not port UNU.RAN's performance-optimized polynomial PINV tables
  one-for-one.
- **AROU/SROU/TDR/ITDR/TABL**: the Fortran constructors and semantics are
  present.  For arbitrary distributions their robust fallback is the same
  numerical inversion engine; thus the target is preserved but the original
  rejection-envelope implementation and performance profile are not copied
  one-for-one.
- **DARI/DAU/DGT**: discrete generator constructors use exact discrete
  inversion (or direct named-distribution RNGs).  The upstream guide/alias
  table optimizations are not copied one-for-one.
- **MIXT**: direct mixture selection followed by component sampling.
- **HITRO**: custom multivariate targets use a hit-and-run-direction symmetric
  Metropolis kernel.  Named multivariate normal/t/Cauchy laws use direct
  generation.
- **VNROU**: named multivariate laws use direct generation; for an arbitrary
  callback target the current v0.1 implementation shares the HITRO
  target-preserving Markov fallback rather than UNU.RAN's independent
  multivariate ratio-of-uniforms rejection kernel.

The distinctions above are deliberate and explicit: this v0.1 package aims
for a useful, fully native Fortran computational API without claiming a
line-for-line port of every internal UNU.RAN optimization.

## Specialized distribution sampling

Several difficult named laws avoid brute-force CDF inversion:

- GIG is transformed with `y = log(x)`.  Its transformed log density is
  strictly concave, and an exact two-tangent rejection envelope is used.
- generalized hyperbolic and hyperbolic variates use the normal-GIG mixture
  representation;
- variance-gamma uses the normal-gamma mixture representation;
- Planck uses the exact zeta mixture of gamma distributions;
- slash uses the normal/uniform ratio representation;
- power-exponential uses a signed gamma-power representation.

Bessel-K evaluation is therefore needed primarily for density evaluation and
reference probability work, not for these samplers.

## R-specific code intentionally omitted

The S4 classes, R parser/string interface, `.Call` registration, packed-object
serialization, R option handling, R RNG-state integration, help/printing code,
and deprecated R aliases are not translated.  The corresponding statistical
objects are explicit Fortran derived types and constructors.

## Validation

The included tests cover standard continuous and discrete reference values,
truncation, custom PDF/CDF/log-PDF callbacks, exact finite probability-vector
distributions, ARS sampling, mixtures, multivariate generation/HITRO,
GIG/GHYP/VG/Meixner/Planck fixed density references, and GIG/GHYP/VG sampling
moments.
