# Porting notes

## Scope

`bridgedist` 0.1.3 has four substantive exported numerical routines:
`dbridge`, `pbridge`, `qbridge`, and `rbridge`. All four are translated.
The R `.onAttach` startup messages and vignette/graphics material are not
computational and are not part of the Fortran API.

## Stable CDF evaluation

Upstream evaluates

`atan((exp(phi*q) + cos(pi*phi))/sin(pi*phi))`.

That can overflow for large positive `q`. The Fortran implementation solves the
same inverse-CDF identity with `atan2` and scales by `exp(-abs(phi*q))`, so it
remains finite far into either tail.

## Stable density evaluation

The direct denominator `cosh(phi*x) + cos(pi*phi)` can overflow. The Fortran
implementation evaluates its logarithm with a scaled expression for large
`abs(phi*x)` and exponentiates only when an ordinary density is requested.

## Tail options

The bridge distribution is symmetric. Therefore `P(X > q) = F(-q)` and the
upper-tail quantile is `-Q(p)`. These identities are used to avoid cancellation.

## Parameter validation

The package documentation requires `0 < phi < 1`, although the R routines do
not explicitly check it. The Fortran routines return IEEE NaN for invalid
`phi` or invalid probabilities rather than propagating arbitrary trigonometric
results.

## R vector recycling

Fortran elemental procedures support scalar expansion but not R's arbitrary
unequal-length recycling. The port therefore provides explicit
`dbridge_recycle`, `pbridge_recycle`, and `qbridge_recycle` helpers.
