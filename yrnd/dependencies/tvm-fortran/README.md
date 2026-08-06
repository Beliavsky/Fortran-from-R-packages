# tvm-fortran

A modern Fortran translation of the computational routines in the R package
`tvm` 0.5.2 by Juan Manuel Truppia.

The library provides time-value-of-money calculations, loan cashflow
construction, internal-rate-of-return solvers, discount/rate transformations,
and interpolated interest-rate curves. It is self-contained and has no external
library dependencies.

## License

`tvm` declares `MIT + file LICENSE`. This translation preserves the MIT license
and the original 2014 copyright attribution. See `LICENSE`, `NOTICE`, and the
unmodified package retained under `original/tvm-0.5.2`.

## Build

```text
fpm build
fpm test
fpm run tvm_demo
fpm run --example basic_cashflows
fpm run --example rate_curves
```

## Main modules

- `tvm`: convenient public umbrella module.
- `tvm_cashflows`: NPV, IRR, loans, payments, balances, and cashflow discounting.
- `tvm_curves`: rate conversions, curve construction, interpolation, and valuation.
- `tvm_interpolation`: monotone piecewise-cubic Hermite interpolation.
- `tvm_root`: self-contained scalar root finding.
- `tvm_kinds`: the `dp = kind(1.0d0)` real kind.

## Cashflow example

```fortran
program cashflow_example
   use tvm
   implicit none
   type(loan_t) :: mortgage
   real(dp) :: value

   mortgage = loan(0.05_dp, 10, 100.0_dp, "french")
   value = npv(0.01_dp, [-1.0_dp, 0.5_dp, 0.9_dp], &
      [0.0_dp, 1.0_dp, 3.0_dp])

   print *, mortgage%cf
   print *, value
   print *, irr([-1.0_dp, 0.5_dp, 0.9_dp], &
      [0.0_dp, 1.0_dp, 3.0_dp])
end program cashflow_example
```

The generic `xnpv` and `xirr` interfaces accept either real year fractions or
integer day serials. Integer dates use ACT/365 fractions relative to the first
date. Use `continuous_compounding` as the compounding-frequency argument for
continuous compounding and `0.0_dp` for simple interest.

## Rate-curve example

```fortran
program curve_example
   use tvm
   implicit none
   type(rate_curve_t) :: curve
   real(dp), allocatable :: swap_rates(:)

   curve = rate_curve_from_rates([0.10_dp, 0.20_dp, 0.30_dp], &
      "zero_eff")
   swap_rates = curve%rate_grid("swap")

   print *, curve%discount(curve%knots)
   print *, swap_rates
   print *, curve%present_value([-1.0_dp, 1.10_dp], &
      [0.0_dp, 1.0_dp])
end program curve_example
```

Curve constructors are available for:

- vectors of rates and maturities;
- vectors of discount factors and maturities;
- a scalar rate callback sampled at supplied knots;
- a scalar discount-factor callback retained for direct evaluation.

Supported rate type strings are `fut`, `zero_nom`, `zero_eff`, `swap`,
`zero_cont`, `french`, and `german`. The last two are output curve types; they
are not direct constructor inputs, matching the original package.

## Deliberate Fortran API choices

- R S3 loans and curves are represented by `loan_t` and `rate_curve_t`.
- R `Date` objects are represented by integer day serials.
- R function objects are represented by explicit Fortran procedure callbacks.
- Plotting is not compiled. Curves return all numerical values needed by any
  plotting library or output program.
- Invalid dimensions and unsupported rate types fail early with descriptive
  messages.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.
