# SmithWilsonYieldCurve modern Fortran

A self-contained modern Fortran translation of the computational code in the
R package `SmithWilsonYieldCurve` 1.1.1.

The library calibrates a Smith-Wilson discount curve to zero-coupon cashflows
or to simple LIBOR, fixed-rate swap, and fixed-rate bond instruments. It uses
only standard Fortran and has no BLAS, LAPACK, R, or C dependencies.

## Implemented functionality

- Smith-Wilson kernel evaluation
- Symmetric Wilson kernel-matrix construction
- Calibration of kernel weights by pivoted Gaussian elimination
- Discount-price and continuously compounded spot-rate evaluation
- Exact repricing of calibration instruments
- Compound-kernel evaluation corresponding to the upstream `K` closure
- LIBOR, swap, and bond payment schedules and cashflows
- Construction of the combined cashflow-time vector and cashflow matrix
- Convenience fitting directly from typed market instruments
- Explicit status codes and diagnostic messages

Plotting and R-specific closures, lists, data frames, S3 classes, and package
hooks are intentionally omitted.

## Build with FPM

```text
fpm test
fpm run --example example_smith_wilson
```

The package name in `fpm.toml` is `smith_wilson_yield_curve`.

## Build with GNU Make

Checked build with runtime bounds and validity checks:

```text
make check
```

Optimized build:

```text
make release
```

Windows command scripts using `gfortran` are supplied in `scripts/`.

## Minimal example

```fortran
program fit_curve
   use smith_wilson_yield_curve, only : dp, fit_smith_wilson_curve_to_instruments, &
                                        make_market_instrument, market_instrument, &
                                        smith_wilson_curve, sw_success
   implicit none

   type(market_instrument) :: instruments(2)
   type(smith_wilson_curve) :: curve
   integer :: info
   character(len=256) :: message

   call make_market_instrument('SWAP', 1.0_dp, 0.025_dp, instruments(1), &
                               info, message, frequency=1.0_dp)
   call make_market_instrument('SWAP', 10.0_dp, 0.05_dp, instruments(2), &
                               info, message, frequency=1.0_dp)

   call fit_smith_wilson_curve_to_instruments(instruments, 0.04_dp, 0.1_dp, &
                                               curve, info, message)
   if (info /= sw_success) then
      print '(a)', trim(message)
      error stop 1
   end if

   print '(f12.8)', curve%discount(20.0_dp)
   print '(f12.8)', curve%continuous_spot(20.0_dp)
end program fit_curve
```

The complete example is in `example/example_smith_wilson.f90`.

## Main API

The facade module is:

```fortran
use smith_wilson_yield_curve
```

Important procedures and types are:

- `type(market_instrument)`
- `type(smith_wilson_curve)`
- `make_market_instrument`
- `fit_smith_wilson_curve`
- `fit_smith_wilson_curve_to_instruments`
- `create_time_vector`
- `create_cashflow_matrix`
- `create_market_value_vector`
- `create_kernel_matrix`
- `fit_kernel_weights`
- `wilson_function`

A fitted curve provides the type-bound methods `discount`,
`continuous_spot`, `compound_kernel`, and `repriced_values`.

## Numerical conventions

The translation preserves the upstream convention that `ufr` is a
continuously compounded rate and therefore uses `exp(-ufr*t)` directly.
LIBOR and swap market values are one per unit notional. A bond uses its
`price` component. Swap cashflows require `tenor*frequency` to be an integer.
Bond schedules are generated backward from maturity, preserving the upstream
handling of a short first coupon period.

## Validation

Five deterministic tests cover kernel reference values, the published
zero-coupon example, exact calibration-instrument repricing, all instrument
schedules, mixed-instrument fitting, long-term convergence, singular systems,
and invalid inputs.

See `docs/VALIDATION.md` for details.

## License

The upstream package declares `GPL-3`. This translation is distributed under
`GPL-3.0-only`. The full license is in `LICENSE` and `license/GPL-3.0.txt`.
The complete unmodified upstream source and original ZIP are retained in
`upstream/` for license and provenance purposes.
