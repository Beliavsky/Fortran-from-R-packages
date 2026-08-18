# Translation notes

## Upstream package

- Package: `contfrac`
- Version: 1.1-12
- Author: Robin K. S. Hankin
- Upstream license declaration: `GPL-2`
- Upstream computational sources: `R/contfrac.R`, `src/contfrac.c`

## R/C to Fortran mapping

| Upstream routine | Fortran routine | Notes |
| --- | --- | --- |
| `CF()` | generic `cf` | Real and complex overloads. |
| `GCF()` | generic `gcf` | Modified Lentz iteration; real and complex overloads. |
| `c_contfrac()` | `gcf` real implementation | Native Fortran arithmetic. |
| `c_contfrac_complex()` | `gcf` complex implementation | Native `complex(dp)` arithmetic replaces split real/imaginary C arrays. |
| `convergents()` | generic `convergents` | Returns allocatable numerator and denominator arrays through subroutine arguments. |
| `gconvergents()` | generic `gconvergents` | Same classical three-term recurrence as upstream. |
| `c_convergents*()` | `gconvergents` implementations | Native real/complex Fortran. |
| `as_cf()` | `as_cf` | Same floor/reciprocal algorithm. |

## Numerical behavior retained

The modified Lentz implementation follows the upstream initialization
`f = C = 1e-30`, `D = 0`, and uses the upstream C literal `2.22044604925031e-16` as the internal
stopping criterion.  This is intentionally slightly smaller than the exact
Fortran `epsilon(1.0_dp)` value and matters for the complex tangent regression
test.  The caller tolerance is used to decide whether an unfinished finite
input is acceptable when `finite=.false.`, just as in the R wrapper.

`CF()` upstream treats the first infinite coefficient as a terminator and marks
the resulting fraction finite.  The Fortran `cf` overloads do the same.  A
nonconverged non-finite call returns IEEE quiet NaN and optionally reports the
residual/iteration count through `contfrac_info`.

The classical convergent recurrence can overflow for long continued fractions,
just as noted in the upstream documentation.  `cf`/`gcf` should be used when
only the limiting value is required.

## Deliberate Fortran API differences

- R lists are replaced by allocatable output arrays for numerators and
  denominators.
- R warnings are replaced by an IEEE NaN result plus optional `contfrac_info`.
- Invalid dimensions are programming errors and use `error stop`.
- `as_cf` returns a length-`n` real array.  For an exactly representable rational
  value, it terminates at the exact final coefficient and fills the remaining
  entries with quiet NaN; optional `n_used` reports the valid prefix length.
  The upstream R function instead continues into division by zero and documents
  that rational values are problematic.

## Omitted code

There is no plotting routine in the package's computational source.  R-specific
`.C` registration, argument coercion, list construction, warning machinery, and
package namespace glue are not translated because they are not numerical
algorithms.
