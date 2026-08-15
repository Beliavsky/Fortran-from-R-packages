# Translation notes

## Scope

`rmutil` contains two different kinds of code: numerical kernels and extensive
R infrastructure for repeated-measurement objects/formulas. Version 0.1.0
ports the numerical kernels and intentionally does not emulate R's S3 or formula
runtime.

## TOMS 614

Upstream ships an f2c-derived C wrapper/source for INTHP, ACM TOMS Algorithm 614.
The Fortran module follows the original algorithm statements and exposes an
explicit procedure interface. The historical algorithm retains its original
control flow where doing so makes provenance and comparison easier.

## Numerical integration changes

The upstream `int2` nests its vectorized Romberg engine. `integrate_2d` instead
uses product Gauss-Legendre quadrature with mappings for semi-infinite and
infinite limits. The requested tolerance selects a 24-, 48-, or 72-point rule
per dimension. This avoids global callback state and remains thread-friendly.

The generalized inverse-Gaussian normalizing Bessel K function is evaluated
from its integral representation, so no external BLAS/LAPACK/special-function
library is required.

## Matrix exponential

Upstream `mexp` offers spectral decomposition or a finite series approximation.
The port supplies a self-contained scaling-and-squaring Taylor matrix
exponential. This is sufficient for `lin.diff.eqn` and avoids an eigensolver
runtime dependency.

## R object layer omitted

`restovec`, `tcctomat`, `tvctomat`, `dftorep`, `rmna`, `lvna`, profile methods,
formula interpreters, print methods, plotting, and readers mostly validate,
reshape, annotate, or dispatch R objects. They are not part of the Fortran
numerical API. Plain arrays plus `nobs`/`nknt` vectors replace those objects in
computational calls such as `gettvc`.

## Compatibility quirks

- Multiplicative binomial: density and CDF paths intentionally follow their
  respective upstream C formulas even though those formulas are not identical.
- PVF Poisson: the upstream CDF's strict `j < q` loop is corrected to the
  documented inclusive cumulative convention. Quantiles use the corrected CDF.
- Skew Laplace: the upstream quantile/RNG branch at `p=0.5` is preserved.
- `mu1.1o2cc`: the formula is translated literally, including the upstream
  final denominator `(beta-ka)^2`.
