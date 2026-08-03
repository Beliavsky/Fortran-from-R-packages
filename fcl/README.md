# fcl-fortran

A modern Fortran translation of the computational core of the R package
[`fcl`](https://github.com/shrektan/fcl), packaged for the Fortran Package
Manager (FPM).

The translation includes:

- calendar-date arithmetic used by the original package;
- dated net present value (`xnpv`) and internal rate of return (`xirr`);
- fixed-rate and zero-coupon bond cash-flow schedules;
- yield to maturity, Macaulay duration, and modified duration;
- daily and cumulative time-weighted returns;
- cumulative profit and loss;
- Modified Dietz average capital and returns.

The R, `extendr`, and `xts` interface layers are not translated. The Fortran API
uses derived types and allocatable arrays directly.

## Build and run

```sh
fpm build
fpm test
fpm run
```

## Main API

```fortran
use fcl

type(fixed_bond_type) :: bond
type(bond_value_type) :: result

bond%value_date = make_date(2021, 2, 1)
bond%maturity_date = make_date(2030, 2, 1)
bond%redemption_value = 100.0_dp
bond%coupon_rate = 0.03_dp
bond%coupon_frequency = 1

result = bond%value(make_date(2022, 2, 1), 100.0_dp)
```

For return calculations, dates are integer day indices, matching the original
R implementation's internal representation of `Date` values.

## Numerical and behavioral notes

- Supported coupon frequencies are 0, 1, 2, 4, 6, and 12.
- Frequency 0 means a single coupon at maturity, matching the source package.
- Coupon accrual uses actual calendar days within each coupon period.
- Discount exponents use the original package's `year_frac` formula, which is
  similar to a 30/360 convention but retains a day adjustment divided by 365.
- Invalid calculations return a nonzero status and, for bond values, IEEE NaNs.
- No business-day calendar or weekend adjustment is applied.

## License

MIT License, preserving the original package license and copyright notice.
See `LICENSE`.
