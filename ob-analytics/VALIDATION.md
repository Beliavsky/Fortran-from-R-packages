# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64 validation environment
- Validation date: 2026-07-25

FPM itself was not installed in the validation environment. The project was
therefore compiled directly in the same source dependency order represented by
the FPM module graph. `fpm.toml` was parsed independently as TOML and audited
for automatic `src`, `app`, `example`, and `test` discovery.

## Strict checked build

The included `scripts/validate.sh` uses:

```text
-std=f2018
-Wall
-Wextra
-Wconversion-extra
-Wimplicit-interface
-Wno-compare-reals
-Werror
-fcheck=all
-fbacktrace
-O0
```

`-Wno-compare-reals` is intentional. Exact comparisons are used for rounded
CSV prices and volumes, fill-size grouping, and consecutive book-state change
detection, matching the original package semantics. All other enabled warnings
are promoted to errors.

Results:

```text
test_alignment: PASS
test_depth_book: PASS
test_order_types: PASS
test_processing: PASS
test_trades: PASS
validation: PASS
```

## Optimized build

All modules, tests, the application, and both examples were rebuilt with `-O2`
and the same static warning checks. Results:

```text
test_alignment: PASS
test_depth_book: PASS
test_order_types: PASS
test_processing: PASS
test_trades: PASS
optimized: PASS
```

## Tested behavior

- Needleman-Wunsch alignment for the documented `2,4,5` versus `1:5` example.
- The original package's simple and conflicting event-matching cases.
- Ambiguous equal-fill bursts and cutoff enforcement.
- Maker/taker assignment from exchange timestamp and order ID.
- Trade direction, maker price, grouped impact volume, hit count, and VWAP.
- Flashed-limit, resting-limit, pacman, market, and market-limit classification.
- Price-level cumulative depth after additions, partial fills, and deletions.
- Best bid/ask and best-level volume through time.
- Filtered depth opening and closing clamps.
- Spread-state deduplication.
- Instantaneous order-book reconstruction and latest changed volume.
- CSV negative-volume removal, life-cycle ordering, fill calculation, matching,
  classification, depth generation, and end-to-end processing.
- Demo and both example executables.

## Structural audits

```text
manifest: PASS
ascii: PASS
source audit: PASS
```

The release contains:

- 11 library source modules
- 5 test programs
- 1 application
- 2 examples
- 19 Fortran files
- 2,253 lines of Fortran

Every translated Fortran file:

- carries `SPDX-License-Identifier: GPL-2.0-or-later`;
- contains `implicit none`;
- is ASCII text;
- has no line longer than 132 columns.

## Runtime defect prevention

Checked builds found three early translation sites that had incorrectly relied
on C/R-style short-circuit evaluation. Fortran does not require short-circuit
logical evaluation. Those expressions were rewritten as nested conditionals,
and the complete suite now passes with `-fcheck=all`.
