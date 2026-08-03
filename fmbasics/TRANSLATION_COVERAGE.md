# Translation coverage

## Included computational areas

| Upstream area | Fortran representation |
|---|---|
| Currency and benchmark currency constructors | `currency_t`; `aud`, `usd`, etc. |
| Currency pairs and inversion | `currency_pair_t`; named pair constructors, `invert` |
| FX spot, today, tomorrow and forward dates | `to_spot`, `to_fx_value`, related routines |
| Ibor and cash indices | `index_t`; all exported benchmark constructors |
| Interest rates and discount factors | Typed constructors, conversions and arithmetic |
| Interpolation objects | Constant, linear, cubic, log-DF and time-variance methods |
| Zero curves | Construction, CSV loading, zeros, discount and forward rates |
| CDS specifications and quotes | Typed CDS specification and curve objects |
| Survival and hazard rates | Conversion, arithmetic and interpolation |
| CDS curve construction | Internal piecewise-hazard bootstrap |
| Volatility quotes and surfaces | CSV loading and two-dimensional interpolation |
| Money and cash flows | Single/multicurrency objects, aggregation and dated flows |

All benchmark constructors exported by the R package are represented, including
AONIA, BBSW, EURIBOR, LIBOR/TIBOR/HIBOR/NIBOR variants, EONIA, SONIA, TONAR,
NZIONA, Fed Funds, CHF TOIS and HONIX.

## Adapted dependencies

- `fmdates`: replaced by internal date, day-count and calendar routines.
- `credule`: replaced by a self-contained CDS survival bootstrap.
- `stats::splinefun`: replaced by an internal natural cubic spline.
- `readr`/`tidyr`: replaced by small readers for retained package CSV data.

## Omitted R-specific infrastructure

- S3 registration and class predicates
- print and format methods
- tibble conversion and summary display
- package-load hooks and R namespace machinery
- documentation-site and vignette generation

These omissions do not remove numerical algorithms. Differences arising from
the adapted calendar database and CDS schedule/bootstrap conventions are
explained in `PORTING.md`.
