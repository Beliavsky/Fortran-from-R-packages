# Porting notes

## Representation

The original package stores the option-value and hedge arrays in C arrays and
then flattens them in R column-major order. The Fortran port stores them directly
as `c(n_periods, n_grid)` and `a(n_periods, n_grid)`, which is the matrix shape
seen by R users.

The original integer `put` flag is represented by the logical argument
`is_put`. `.true.` selects a put and `.false.` selects a call.

## Discounted quantities

The algorithm works with discounted asset prices and uses the discounted strike

```text
strike * exp(-rate * maturity)
```

This behavior is unchanged. Input returns are periodic log excess returns. The
port computes `xi = exp(log_excess_returns) - 1` exactly as the C routine did.

## Interpolation

The original routine linearly interpolates on a uniform grid and linearly
extrapolates beyond both endpoints. The Fortran implementation preserves that
behavior. It additionally handles a one-point vector safely.

## Safety corrections

The numerical formulas are unchanged, but several implementation defects and
undefined cases in the C source are corrected:

1. The C pointer arrays were allocated with `sizeof(double)` instead of
   `sizeof(double *)`. This happens to work on many 64-bit systems but is not
   portable. Fortran allocatable arrays remove the issue.
2. The temporary `s0` allocation was never freed. The Fortran implementation
   has no corresponding leak.
3. The C routine did not validate sample size, grid size, grid bounds, the
   return second moment, or the singular case `1 - E[xi]^2/E[xi^2] = 0`.
   The Fortran result reports these failures through `ok` and `message`.
4. The original initial-share calculation divides by every grid price without
   checking for zero. The translated hedging engine requires a strictly
   positive price grid.
5. The R/C interface performed several avoidable scalar heap allocations.
   The Fortran implementation uses scalar local variables.

These changes affect invalid or nonportable cases; valid-input numerical output
is preserved.

## Random numbers

The original package relied on R to simulate returns. The core function still
accepts caller-supplied returns. A small native Box-Muller normal generator is
included only to make examples self-contained and reproducible.

## Licensing

The package DESCRIPTION declares `GPL (>= 2)`. The FPM manifest and all new
Fortran files use the SPDX identifier `GPL-2.0-or-later`. The complete GPLv2
text is included, and the original package is retained unmodified.
