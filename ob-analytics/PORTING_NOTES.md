# Porting notes

## Language and storage choices

- Timestamps are `integer(int64)` milliseconds instead of POSIXct objects.
- Prices and volumes use `real(real64)`.
- R factors are integer constants.
- Data frames and lists are typed derived structures with allocatable arrays.
- The implementation is self-contained and uses no BLAS, LAPACK, or statistics
  dependency.

## Sparse depth state

The R `depthMetrics` implementation allocates two arrays of one million price
slots and indexes them by integer cents. The Fortran port keeps sparse arrays
of active integer-cent levels and volumes. This removes the hard-coded price
ceiling while preserving basis-point bin calculations.

## Sequence matching

`event_match` preserves the original two-stage method:

1. For each common fill size, match each bid fill to its nearest ask fill.
2. If candidate matches collide, use Needleman-Wunsch global alignment with
   the original time-distance score and gap penalty.

Fortran does not guarantee short-circuit evaluation. All potentially unsafe
indexing conditions are therefore written as nested conditionals. This was
explicitly exercised with runtime bounds checking.

## Large price-jump correction

The upstream `matchTrades` correction indexes `diff(price)` locations in a way
that can refer to the record before the actual jump. The Fortran implementation
checks each current trade against its immediately preceding trade and swaps the
current maker/taker roles when the threshold is exceeded. The correction can be
disabled with `correct_large_jumps=.false.` and the threshold is configurable.

## Empty and degenerate values

- `normalize` returns zeros for a constant vector instead of dividing by zero.
- VWAP functions return zero when total weight is zero.
- Depth and book routines return allocated zero-length arrays when no records
  satisfy a filter.
- Missing event matches are represented by integer zero rather than R `NA`.
- Missing aggressiveness values use `has_aggressiveness=.false.`.

## CSV scope

The parser intentionally targets the documented seven-column package format.
It is not a general quoted-field CSV parser. Actions named `modified` are
accepted as aliases for the R implementation's `changed` level.

## Licensing

The package DESCRIPTION declares `GPL (>= 2)`, so the translation uses the SPDX
identifier `GPL-2.0-or-later`. The original package and full GPLv2 text are
included unchanged.
