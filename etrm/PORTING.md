# Porting notes

## Source basis

This project translates the computational content of etrm 1.0.2. The source
package declares `MIT + file LICENSE`; the resulting Fortran source is also MIT
licensed and carries SPDX headers plus the original copyright attribution.

## R-to-Fortran mapping

| R entry point | Fortran entry point | Status |
| --- | --- | --- |
| `cppi` | `cppi` | translated |
| `dppi` | `dppi` | translated |
| `obpi` | `obpi` | translated |
| `shpi` | `shpi` | translated |
| `slpi` | `slpi` | translated |
| `msfc` | `msfc` / `maximum_smoothness_forward_curve` | translated |
| S4 strategy classes | `strategy_result` | replaced by typed value |
| S4 MSFC class | `msfc_result` | replaced by typed value |
| strategy `summary` | `summarize_strategy` | translated numerically |
| plots and data-frame display | none | omitted |

## Numerical representation

- R `Date` values are represented by integer day offsets relative to the trade
  date. The curve time scale remains days divided by 365, exactly as in the R
  routine.
- The MSFC KKT system is solved by LAPACK `dgesv` rather than R `solve`.
- Quartic coefficients use the same order as the source: `a,b,c,d,e`.
- Daily MSFC values retain the source's average over a `0.0001`-day interval,
  rather than evaluating only at the left endpoint.
- Integer trade restrictions use round-to-nearest, ties-to-even behavior to
  mirror R's `round(x, 0)`. Unrestricted trades retain full double precision
  instead of being rounded to ten decimal places.

## Defensive corrections

The following edge cases are handled more explicitly than in the R source:

1. Zero hedge volume is rejected before divisions by `q`.
2. CPPI/DPPI risk percentages must be positive.
3. OBPI and SHPI require `days_left` to cover every supplied observation. The R
   code otherwise generates unavailable times or hedge fractions beyond expiry.
4. SLPI and SHPI hedge immediately when the cap/floor is first reached. In the
   R implementation, a crossing at the first observation updates the hedge
   fraction but not the first position because recomputation begins at index 2.
5. MSFC contracts cannot start before the trade date, and a zero-day knot is
   stored only once. This prevents negative or zero-width spline intervals.
6. The unused R expression `a <- setdiff(S, x)` is omitted. `setdiff` is not a
   valid way to extract Lagrange multipliers when numeric values repeat, but the
   resulting object is never used.
7. The ineffective R length test `length(Date > length(MSFC))` and the unused
   `CompAvg` calculation are omitted.

## Scope exclusions

Plotting, S4 dispatch, data frames, bundled synthetic datasets, and vignette
rendering are R infrastructure rather than numerical algorithms. Original R
source and metadata remain under `original/`, but binary `.rda` datasets are not
redistributed in the Fortran project.
