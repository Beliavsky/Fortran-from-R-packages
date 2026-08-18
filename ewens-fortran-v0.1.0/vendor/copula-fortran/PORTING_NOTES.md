# Porting notes

## Typed models instead of S4 objects

R constructors return S4 instances with slots and generic method dispatch. The
Fortran constructors return a `copula_model` value containing a family code,
dimension, scalar parameters, rotation, and an optional correlation matrix.
Invalid models can be checked with `model%valid()`.

## Probability calculations

Bivariate Gaussian probabilities are calculated by deterministic Simpson
integration of the conditional-normal identity. Bivariate Student probabilities
use the corresponding conditional Student identity.

For dimension greater than two, Gaussian probabilities use a conditional
shifted-Halton integral. Student probabilities use chi-square mixing around the
same Gaussian engine. These are dependency-free alternatives and are not
bit-for-bit replacements for the upstream `mvtnorm`/Genz-Bretz stack.

## Densities

Gaussian and Student copula densities use closed-form multivariate formulas.
For other absolutely continuous bivariate families, the density is the mixed
second derivative of the CDF evaluated by a central difference. This gives one
consistent implementation for many families while keeping the source compact.

Marshall-Olkin and Frechet bounds contain singular mass. Their ordinary density
is not a complete representation of the distribution, so `dCopula` returns zero
for those singular implementations.

## Random generation

- Gaussian and Student copulas use Cholesky factors.
- Clayton uses gamma frailty.
- Gumbel uses positive-stable frailty.
- Positive-parameter multidimensional Frank uses logarithmic-series frailty.
- Other bivariate families use numerical inversion of the conditional CDF.
- Rotations are applied after base-family simulation.

All routines can use deterministic 64-bit seeds and do not depend on the
processor's intrinsic random-number stream.

## Dependence measures

Closed forms are used for the principal elliptical and Archimedean families.
Other smooth bivariate families use deterministic numerical integration.
Parameter inversion uses analytical formulas where available and bisection
otherwise.

## Estimation

`fit_copula` currently supports one free dependence parameter. The `itau`
method inverts the sample Kendall tau. The `mpl` method maximizes the
pseudo-log-likelihood by bounded golden-section search. Standard errors use a
finite-difference observed information calculation.

This is intentionally smaller than upstream `fitCopula`, which supports partly
fixed parameter vectors, multiple optimizers, rank-based corrections, and many
variance estimators.

## Empirical tests

The included independence, exchangeability, and radial-symmetry tests use
simple Cramer-von-Mises-style empirical-copula statistics with permutation or
randomization calibration. They are useful standalone tests but are not
presented as numerical reproductions of every upstream test statistic.

## Numerical safeguards

- Inputs are constrained to valid copula domains.
- Correlation matrices are checked by Cholesky decomposition before density or
  simulation calculations.
- Probability outputs are clipped to valid bounds to suppress roundoff leakage.
- Central-difference steps shrink near unit-square boundaries.
- Rank-correlation denominators are converted to double precision before
  multiplication, avoiding 32-bit integer overflow for moderate samples.
- Real-valued ties use a tight machine-precision comparison rather than exact
  floating-point equality.
- No logical expression relies on short-circuit evaluation.

## Thread safety

The random generator maintains module-level state. Independent seeded calls are
reproducible, but simultaneous calls from multiple threads should use external
serialization or separate processes.
