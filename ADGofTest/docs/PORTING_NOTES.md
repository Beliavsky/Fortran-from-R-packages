# Porting notes

## Exact computational translation

The statistic and Marsaglia finite-sample correction polynomials are translated
directly from the three R source files. Binary64 constants retain the decimal
values used by the upstream code.

## CDF callbacks

R passes an arbitrary function and `...` arguments. Modern Fortran instead uses
a scalar procedure callback. Parameters can be held in module variables or
captured by a dedicated wrapper procedure.

## Endpoint behavior

The upstream range check accepts 0 and 1. Its logarithm then produces an
infinite statistic, and the final p-value is 0. This remains the default.
`clip_probabilities=.true.` replaces endpoints by the nearest safe binary64
values before evaluating logarithms.

## Approximation range

For very small statistics and finite `n`, the published polynomial correction
can return a CDF slightly below 0, making the source-compatible p-value slightly
above 1. The default preserves this behavior. Set `clamp_p_value=.true.` or call
`ad_distribution_cdf(..., clamp_probability=.true.)` to constrain results to
valid probability bounds.

## Omitted R infrastructure

The port omits `htest` class creation, captured expression names, package
namespace machinery, and R's dynamic ellipsis dispatch. No plotting code exists
in the upstream package.
