# ob_analytics for modern Fortran

`ob_analytics` is a self-contained modern Fortran/FPM translation of the
computational parts of the R package `obAnalytics` 0.1.2. It processes limit
order events, infers executions, reconstructs order-book depth, and calculates
market-microstructure summaries.

The original package is Copyright (C) 2015,2016 Philip Stubbings and is
licensed under GPL version 2 or later. The same license is preserved here.

## Implemented functionality

- CSV event loading, rounding, negative-volume removal, duplicate removal, and
  order-life-cycle sorting.
- Per-event fill calculations.
- Bid/ask fill matching using nearest timestamps and Needleman-Wunsch sequence
  alignment for ambiguous bursts.
- Maker/taker inference and trade reconstruction.
- Flashed-limit, resting-limit, market-limit, pacman, market, and unknown order
  classification.
- Price-level cumulative depth through time.
- Time-window depth filtering with clamped opening and closing states.
- Best bid/ask, best-level volume, and configurable basis-point depth bins.
- Spread-change extraction.
- Order aggressiveness relative to the preceding best quote.
- Instantaneous order-book reconstruction with price-level, BPS, liquidity,
  depth-range, and maximum-level filters.
- Market-impact summaries grouped by taker order.
- End-to-end `process_data` pipeline and zombie-order filtering.
- Utility routines for interval sums, VWAP, normalization, reversed matrices,
  and price-level gaps.

Plotting, `zoo` objects, `ggplot2` themes, and RDS serialization are not
compiled. All arrays needed for plotting are returned to the caller. The
original R source and bundled data remain in `original/`.

## Build and test

```text
fpm build
fpm test
fpm run ob_analytics_demo
fpm run --example event_matching
fpm run --example depth_and_book
```

The library has no external dependencies.

## Data representation

Timestamps are signed 64-bit integer milliseconds. Prices and volumes use
`real(dp)`, where `dp = real64`. Categorical R factors are integer constants:

```fortran
use ob_analytics

integer :: action, side, order_type
action = action_created
side = side_bid
order_type = type_resting_limit
```

The central records are:

- `event_t`
- `trade_t`
- `depth_update_t`
- `depth_summary_t`
- `spread_t`
- `order_level_t`
- `order_book_t`
- `impact_t`
- `processing_result_t`

## End-to-end CSV processing

The input CSV has the seven columns expected by the original package:

```text
id,timestamp,exchange.timestamp,price,volume,action,direction
```

Timestamps are integer milliseconds. Actions may be `created`, `changed`,
`modified`, or `deleted`; `modified` is treated as `changed`.

```fortran
use ob_analytics
implicit none

type(processing_result_t) :: result
integer :: status
character(len=:), allocatable :: message

call process_data('orders.csv', result, status=status, message=message)
if (status /= 0) error stop message

print *, size(result%events)
print *, size(result%trades)
print *, size(result%depth)
```

The default pipeline removes inferred zombie orders and removes the first
60,000 milliseconds of depth summaries, matching the original warm-up rule.
Both behaviors are configurable.

## Direct API examples

```fortran
call event_match(events, cutoff_ms=5000_i8)
trades = match_trades(events)
call set_order_types(events, trades)
depth = price_level_volume(events)
summary = depth_metrics(depth, bps=25, bins=20)
spread = get_spread(summary)
book = reconstruct_order_book(events, timestamp_ms)
impacts = trade_impacts(trades)
```

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for detailed mapping,
behavioral differences, and test information.
