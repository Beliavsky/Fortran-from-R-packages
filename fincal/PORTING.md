# Porting notes

## Preserved behavior

The formulas and sign conventions follow FinCal 0.6.3. The translation covers all 54 computational functions exported by the R package.

R names containing periods are represented with underscores. The original misspelling `r.norminal` is represented both by the correctly named `nominal_rate` and by compatibility alias `r_norminal`.

## Root finding

R's `uniroot` is not available in a self-contained Fortran package. `discount_rate` and `irr` therefore use:

1. a scan equally spaced in `log(1 + rate)`, which gives useful resolution near zero while still covering large positive rates; and
2. bisection after a sign-changing interval is found.

`irr` searches nonnegative rates by default, matching the practical role of the original `irr`. `solve_irr` accepts explicit lower and upper limits. `irr2` preserves the original threshold-based grid search for negative or otherwise difficult roots.

Projects with non-conventional cash flows can have multiple IRRs. These routines return the first root found in increasing rate order within the requested interval.

## Zero-rate limits

The R formulas for annuities and payments divide by `r`. The Fortran implementation evaluates the mathematically correct limits when `r` is numerically zero.

## Array and error handling

R exceptions are replaced by status codes and IEEE quiet NaNs where a scalar result cannot be returned. Inventory costing uses `type(inventory_result)`. Double-declining depreciation returns an allocatable vector containing one expense per life year.

## Inventory costing

FIFO, LIFO, and weighted-average costing preserve the original definitions. The Fortran implementation uses a common lot-depletion algorithm, avoiding the duplicated index logic in the R source while producing the same results for valid inputs.

## Omitted code

Plotting functions and network-based OHLC download wrappers were not translated. They do not contain portable numerical algorithms and would require graphics or HTTP dependencies contrary to the self-contained FPM design.
