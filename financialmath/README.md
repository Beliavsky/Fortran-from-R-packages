# financialmath-fortran

A modern Fortran translation of the computational core of the R package
`FinancialMath` 0.1.1 by Kameron Penn and Jack Schmidt.

The project is an FPM package and has no external numerical dependencies.
Plotting and R-specific matrix/list presentation are deliberately excluded;
typed derived types and arrays return the numerical results needed by callers.

## Build

```text
fpm build
fpm test
fpm run financialmath_demo
fpm run --example cashflows_and_loans
fpm run --example options_and_forwards
```

A direct GNU Fortran validation script is also supplied:

```text
./scripts/validate.sh
```

On Windows, with FPM on `PATH`:

```text
scripts\validate.bat
```

## Main modules

- `financialmath_cashflows`: NPV, IRR, TVM, rate conversion, cash-flow
  duration/convexity, swaps, and dollar/time-weighted yields.
- `financialmath_annuities`: level, arithmetic, and geometric annuities and
  perpetuities, including solving for common unknown parameters.
- `financialmath_loans`: amortization schedules, period decomposition, and
  bond price/duration/convexity calculations.
- `financialmath_derivatives`: Black-Scholes values and first-order Greeks,
  forwards, prepaid forwards, and option-strategy payoff/profit tables.
- `financialmath`: umbrella module re-exporting the public API.

All real calculations use:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Example

```fortran
use financialmath
implicit none

type(tvm_result_t) :: tv
type(option_order1_t) :: bs

tv = solve_tvm(1000.0_dp, 0.0_dp, 10.0_dp, 0.05_dp, 1.0_dp, 'fv')
bs = bls_order1(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp, 0.0_dp)

print *, tv%future_value
print *, bs%call_price, bs%call_delta
```

## API conventions

R functions that represented an unknown with `NA` accept an explicit `unknown`
string in the Fortran API. For example, `annuity_level(..., 'payment')` solves
for the payment and `solve_tvm(..., 'rate')` solves for the interest rate.

The annuity argument `immediate=.true.` means payments occur at the ends of
periods, matching the upstream package's `imm=TRUE`; `.false.` represents an
annuity-due.

For forward contracts, `force_rate` is continuously compounded, matching the
upstream use of `exp(r*t)`. A negative `growth_rate` such as `-1.0_dp` selects
level discrete dividends.

## License and provenance

The upstream package declares `GPL-2`. This translation preserves it as
`GPL-2.0-only`. The complete GPL version 2 text is in `LICENSE`.

The supplied R source is retained unmodified in
`original/FinancialMath-0.1.1`, and the supplied ZIP is retained in
`provenance/FinancialMath-master.zip`. Checksum manifests are included for the
source archive, original files, and translated release files.
