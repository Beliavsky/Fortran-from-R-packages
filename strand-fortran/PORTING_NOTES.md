# Porting notes

## Design

The R package represents cross sections and simulation detail as data frames and
coordinates work through mutable R6 objects. The Fortran port instead exposes
purely numerical arrays, explicit configuration structures, and returned result
structures. This makes dimensions, ownership, and missing-value conventions
visible to the caller and removes dependencies on Rglpk, Matrix, Arrow, and the
tidyverse.

All real calculations use `real(dp)`, where `dp = kind(1.0d0)`. Source is
Fortran 2018, uses `implicit none`, and stays within the standard 132-column
free-form limit.

## Linear programming

`PortOpt` delegates to Rglpk and can optionally use Rsymphony. The port uses a
self-contained bounded two-phase simplex implementation. Signed trade variables
are shifted by their lower bounds, finite upper bounds become constraints, and
phase I supplies artificial variables where necessary.

The objective and constraint construction follows the upstream model:

- one signed dollar trade per security and strategy;
- one nonnegative auxiliary variable per security for absolute joint trade;
- exact target long and short market values;
- position and trade bounds based on expected dollar volume and strategy size;
- factor/category bands, turnover, and optional progressive loosening.

The optimizer reports failure instead of raising an R condition. Inspect
`optimization_result%success` and `%message`.

## Share rounding

The LP is solved in dollar space. As upstream does, each strategy trade is then
rounded to whole shares and its dollar value is recomputed. Therefore a rounded
portfolio can differ slightly from the continuous LP constraints, especially
for high-priced securities.

## Internal transfers and market fills

Orders that offset across strategies are transferred internally. Only the
remaining joint order is sent to the market. When a joint order is volume
limited, the market portion is allocated proportionally among same-direction
strategy orders; the residual transfer portion fills completely.

This preserves the upstream distinction between:

- internal shares, created by cross-strategy transfers; and
- external shares, acquired or sold in the market.

## Corporate actions and delistings

A supplied adjustment ratio divides both internal and external shares before the
day's optimization. For example, a ratio of `0.5` doubles shares for a two-for-one
split, matching the upstream convention.

On a delisting day the position is closed, ordinary mark-to-market P&L and costs
are replaced by starting NMV times the supplied delisting return, and ending
shares are zero. `simulate_portfolio` accepts day-by-security delisting arrays.

## Missing values

IEEE quiet NaN represents missing real data. `rank_normal`, cross-section
statistics, and neutralization use complete finite observations. Integer IDs,
category codes, logical investability flags, and shares do not use sentinel
missing values.

## Normalization and neutralization

`rank_normal` preserves the upstream 11-decimal rounding used before rank tie
calculation. `adjust_numeric` reproduces repeated tail-weighted and ordinary
least-squares residualization followed by rank normalization.

The upstream `adjust` accepts factors and categories through formula/data-frame
machinery. The Fortran routine accepts a numeric design matrix. Encode a
categorical variable as dummy columns when categorical neutralization is needed.

## Category exposure layout

`calculate_exposures` stores category exposure in
`category(level_index, category_variable, strategy)`. The original integer code
for each compact level index is returned in `category_level`; the number of
valid levels for category variable `k` is `category_count(k)`. This compact
mapping avoids allocating large gaps for sparse category codes.

## Performance conventions

Returns are daily P&L divided by daily GMV. Annualization defaults to 252 periods.
Drawdown is calculated from the additive cumulative-return series, matching the
upstream reporting convention. Holding period is derived from average turnover
and average GMV using the upstream annualized two-sided-turnover formula.

## Intentional safety changes

- Strategy capital, prices, volumes, dimensions, and adjustment ratios are
  validated before optimization or simulation.
- Fortran control flow does not assume short-circuit evaluation.
- Multi-day simulation uses distinct day and strategy loop variables.
- Optional delisting and financing arrays are dimension checked.
- The native solver returns explicit feasible/optimal/unbounded status.

## License

The upstream package declares `GPL-3`. This translation uses the SPDX identifier
`GPL-3.0-only` and includes the complete GPLv3 text. The unmodified upstream tree
and supplied archive are retained for provenance.
