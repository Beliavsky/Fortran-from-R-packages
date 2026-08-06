# treasuryTR-fortran

A modern Fortran 2018 translation of the computational routines in the R
package `treasuryTR` 0.1.6.

The library converts constant-maturity yield observations into approximate
bond total returns using modified duration, convexity, and yield income. It is
self-contained and has no BLAS, LAPACK, network, or third-party runtime
dependency.

## Implemented routines

- `mod_duration`: modified duration for scalar or conformable array arguments
- `convexity`: bond convexity for scalar or conformable array arguments
- `period_total_return`: one-period return from current and lagged yields
- `total_return`: vector and matrix time-series interfaces
- `carry_forward`: last-observation-carried-forward missing-value treatment
- `percent_to_decimal`: percentage-rate conversion
- `prepare_yields`: preprocessing equivalent to the non-network part of
  upstream `get_yields`

The first observation returned by `total_return` is IEEE NaN because no lagged
yield exists. Matrix input is interpreted as time by series, and each column is
processed independently.

## Source-compatible and stable modes

The upstream formulas are used by default. Pass
`source_compatible=.false.` to use numerically stable `log1p`/`expm1`
equivalents and the analytic zero-yield limits:

- modified duration: `maturity`
- convexity: `maturity * (2*maturity + 1) / 2`

The source-compatible formulas return NaN at exactly zero yield because the R
expressions contain `0/0` cancellation.

## Example

```fortran
program example
  use treasurytr, only : dp, prepare_yields, total_return
  implicit none

  real(dp) :: yields(4), returns(4)

  yields = prepare_yields([4.00_dp, 4.05_dp, 4.03_dp, 4.10_dp])
  returns = total_return(yields, maturity=10.0_dp, scale=12.0_dp, &
    source_compatible=.false.)
end program example
```

## Build

With GNU Make and GNU Fortran:

```text
make test-checked
make test-optimized
make example
```

With FPM:

```text
fpm test
fpm run
```

The checked build uses strict warnings, bounds checking, and uninitialized-real
initialization. The optimized build uses the same language and warning checks
at `-O3`.

## Layout

- `src/`: library modules
- `app/`: demonstration program
- `test/`: deterministic tests
- `scripts/`: Linux/macOS and Windows build helpers
- `provenance/`: complete upstream source and original input archive

## License

MIT. See `LICENSE` and `NOTICE`.
