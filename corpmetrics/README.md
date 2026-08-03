# corpmetrics-fortran

A modern Fortran/FPM translation of the computational algorithms in the R
package **corpmetrics 1.0**.

The library provides typed numerical results for:

- balance-sheet metrics (`balsh`)
- Capital Asset Pricing Model analysis (`capm`)
- zero-growth, Gordon, and differential-growth dividend models (`ddm`)
- fixed-income price and duration (`fis`)
- net present value and internal rate of return (`idm`)
- income-statement metrics (`insta`)
- loan amortization schedules (`loan`)

The original R package returns formatted data frames. This port returns raw
`real(dp)` values in derived types so results can be reused in other numerical
programs. `round2_r` is provided for the two-decimal, ties-to-even rounding used
by the original loan table.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example investment_example
fpm run --example loan_example
```

## Build without FPM

```text
./scripts/test_gfortran.sh strict
./scripts/test_gfortran.sh release
```

The scripts require GNU Fortran and build all tests, applications, and examples
from a clean directory.

## Minimal example

```fortran
program example
    use corpmetrics, only : dp, investment_result, idm, cm_success
    implicit none

    type(investment_result) :: result
    integer :: status

    call idm([-100.0_dp, 120.0_dp], [0.0_dp, 0.10_dp], result, status)
    if (status /= cm_success) error stop 'idm failed'

    print '(a,f10.4)', 'NPV: ', result%npv
    print '(a,f10.4,a)', 'IRR: ', result%irr_percent, '%'
end program example
```

## Dependencies

The library is self-contained and uses only standard Fortran intrinsics.

## License and provenance

The original package declares `GPL (>= 2)`. This translation is therefore
licensed as **GPL-2.0-or-later**. The complete original package tree is retained
under `original/corpmetrics-master/` for provenance.
