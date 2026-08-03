# Metrics modern Fortran

A dependency-free modern Fortran/FPM translation of the computational code in
the R package Metrics 0.1.4.

## Build and run

```text
fpm build
fpm test
fpm run
fpm run --example regression_metrics
```

## Example

```fortran
program example
    use metrics, only : dp, rmse, mae
    implicit none

    real(dp) :: actual(4), predicted(4)

    actual = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    predicted = [1.1_dp, 1.8_dp, 3.2_dp, 3.7_dp]

    print *, rmse(actual, predicted)
    print *, mae(actual, predicted)
end program example
```

The public umbrella module is `metrics`. See `API_MAP.md` and
`PORTING_NOTES.md` for complete coverage and R-to-Fortran differences.
