# Porting notes

## Array orientation

Time is the first dimension and assets are the second dimension throughout:

```text
prices(time, asset)
positions(time, asset)
weights(time, asset)
```

This matches the usual R matrix layout conceptually, although Fortran stores
arrays in column-major order.

## Timestamps and identifiers

Timestamps are `real(dp)` values in `journal_type`; they may represent integer
observation numbers, serial dates, or elapsed time. Instruments and accounts
are integer IDs. Applications can maintain separate string-label tables.

## Missing values

The computational routines expect finite inputs unless a routine explicitly
states otherwise. R's NA propagation and warning conventions are not emulated.

## Portfolio returns

`portfolio_returns_weights` follows the upstream normalized-wealth algorithm:
wealth starts at one, target weights determine holdings at rebalance rows, and
the residual is cash. Contributions are computed from beginning-of-period
holdings.

`portfolio_returns_positions` forward-fills positions between rebalance rows
and computes contributions relative to prior marked portfolio value.

## P/L

`pl_summary` reproduces the upstream closed-trade cashflow identity and average
buy/sell prices. `pl_path` uses the same average-cost logic as `.pl_stats`,
including reversals through zero.

## Unit prices

Rather than matching timestamp classes, `unit_prices` accepts an integer index
for each cashflow. Multiple cashflows at one index share the same computed unit
price and each receives its own issued-unit amount.

## Backtesting

The R implementation permits callbacks that inspect current open/high/low/close
prices, wealth, cash, positions, timestamps, and mutable globals. Such dynamic
callbacks have no direct static Fortran equivalent. The port therefore uses a
target matrix produced before the accounting pass. Strategies needing dynamic
state can generate that matrix in a separate Fortran routine or call the
accounting loop repeatedly.

## Licensing

Because upstream PMwR is GPL-3, this translation is GPL-3.0-only. The upstream
source snapshot is included unchanged under `upstream/`.
