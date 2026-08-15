# Translation notes

## Upstream

- R package: `expint`
- Upstream version: 0.2-1
- Date: 2026-01-26
- Maintainer/primary author: Vincent Goulet
- Numerical kernels: derived in upstream from GNU Scientific Library code by
  Gerard Jungman and Brian Gough

## Translated computational surface

| Upstream R/C | Fortran |
|---|---|
| `expint()` | generic `expint()` / `expint_value` |
| `expint_E1()` / `expint_E1` | `expint_e1` |
| `expint_E2()` / `expint_E2` | `expint_e2` |
| `expint_En()` / `expint_En` | `expint_en` |
| `expint_Ei()` | `expint_ei` |
| `gammainc()` / `gamma_inc` | `gammainc` / `gamma_inc` |
| R recycling | `expint_recycle`, `gammainc_recycle` |

## Exponential integrals

The six upstream Chebyshev tables (`AE11`, `AE12`, `E11`, `E12`, `AE13`,
`AE14`) are retained coefficient-for-coefficient.  Clenshaw evaluation and the
same argument intervals are used.  The large-argument `E2` asymptotic
polynomial and `En`/incomplete-gamma relation are also retained.

## Incomplete gamma

For `a < 0`, the upstream algorithm is translated directly:

- GSL-derived unconditionally convergent continued fraction when `x > 0.25`;
- special recursion for `-0.5 < a < 0` at small `x`;
- downward recurrence from the fractional part of `a` otherwise.

For `a > 0`, upstream calls R's `gamma()` and upper-tail `pgamma()`.  To make
this package self-contained, the Fortran translation instead computes the
upper incomplete gamma with a standard convergent power series for the lower
regularized gamma when `x < a+1`, and a continued fraction for the upper tail
when `x >= a+1`.

## R-specific code omitted

The following are interface rather than numerical-kernel features and are not
ported:

- `.External` dispatch;
- R SEXPs and attribute propagation;
- R NA distinct from IEEE NaN;
- localization/messages/warnings;
- package registration in `init.c`.

The explicit recycling helpers preserve the useful numerical part of R's
vector recycling without introducing R runtime dependencies.
