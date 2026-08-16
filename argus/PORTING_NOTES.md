# Porting notes

## Scope

The computational API of upstream `argus` 0.1.1 is translated:

- `dargus`: density and log-density
- `pargus`: lower/upper CDF and log-probabilities
- `qargus`: lower/upper quantiles and log-probability input
- `rargus`: random generation by inversion or ratio-of-uniforms
- scalar-parameter and varying-parameter random generation
- R-style numeric argument recycling helpers

The upstream test script contains histogram/line plotting. Plotting was not
ported.

## Runuran removal

Upstream inversion generation uses Runuran `pinv.new()` objects and spline
tables to approximate the inverse CDF very quickly. The Fortran port removes
that dependency. It uses the identity

    F(x) = 1 - P(3/2, chi^2 (1-x^2)/2) / P(3/2, chi^2/2)

and numerically inverts the shape-3/2 regularized incomplete-gamma CDF with a
bracketed Newton method. Thus the Fortran inversion path targets the exact
CDF to floating-point solver tolerance instead of reproducing Runuran's
precomputed approximation tables.

## Shape-3/2 incomplete gamma

No external special-function library is needed. For small arguments the code
uses the lower incomplete-gamma series. For larger arguments it uses

    Q(3/2,z) = erfc(sqrt(z)) + 2 sqrt(z) exp(-z) / sqrt(pi)

which is exact for shape 3/2. Log-probability helpers avoid cancellation for
small `chi` and near-tail CDF calculations.

## Ratio of uniforms

The ratio-of-uniforms generator is translated from `src/argus.c`. One
algebraic rewrite is used in the `chi <= 1` setup:

    xm = chi^2 / (chi^2 + 5 + sqrt(chi^2(chi^2+6)+25))

This is algebraically identical to the upstream expression but avoids
subtracting nearly equal floating-point numbers when `chi` is very small.

## Support handling

The Fortran API explicitly enforces the statistical support `[0,1]`:

- density is zero outside the support
- lower CDF is zero below the support and one above it
- upper CDF behaves conversely

Invalid `chi <= 0` arguments to the scalar density/CDF/quantile functions
return IEEE NaN. Random generators reject invalid `chi` with `error stop`.

## Vector semantics

`dargus`, `pargus`, and `qargus` are elemental, so a scalar parameter can be
used with a conformable array directly. The additional `*_recycle`
subroutines reproduce R's recycling semantics when argument arrays have
different lengths.

`rargus_varying` cycles through its `chi(:)` vector if the output sample is
longer than the parameter vector, matching the useful computational behavior
of the R implementation.
