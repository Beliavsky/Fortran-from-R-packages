# vrtest-fortran

`vrtest-fortran` is a self-contained modern Fortran and FPM translation of
Jae H. Kim's R package `vrtest` 1.2. It implements variance-ratio tests and
other tests of the martingale-difference hypothesis without requiring R,
BLAS, LAPACK, or another external library.

## Build

A Fortran 2018 compiler and FPM are required.

```text
fpm build
fpm test
fpm run vrtest_demo
fpm run --example basic_tests
fpm run --example bootstrap_tests
```

Windows users can run `build.bat`; POSIX users can run `./build.sh`.

## Main modules

- `vrtest`: umbrella module for normal use.
- `vrtest_variance_ratio`: variance-ratio, rank/sign, bootstrap, panel, and
  subsampling procedures.
- `vrtest_spectral`: portmanteau, spectral, Dominguez-Lobato, and Chen-Deo
  procedures.
- `vrtest_types`: result derived types.
- `vrtest_utils`: distribution, quantile, rank, random-number, and linear
  algebra support.

All real calculations use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Minimal example

```fortran
program example
   use vrtest
   implicit none
   real(dp) :: returns(8)
   integer :: periods(3)
   type(lmcd_result) :: result
   integer :: i

   returns = [0.01_dp,-0.004_dp,0.007_dp,0.002_dp, &
      -0.009_dp,0.006_dp,0.012_dp,-0.005_dp]
   periods = [2,3,4]
   result = lo_mackinlay(returns,periods)

   do i = 1, size(periods)
      print '(i4,2f14.6)',periods(i),result%homoskedastic(i), &
         result%heteroskedastic(i)
   end do
end program example
```

## Public procedures

### Variance-ratio family

- `ar1_fit`
- `adjust_thin`
- `abel_bandwidth`
- `quadratic_spectral_kernel`
- `fast_variance_ratio`
- `lm_statistic`
- `lmcd_statistics`
- `lo_mackinlay`
- `chow_denning`
- `wald_test`
- `automatic_variance_ratio`
- `variance_ratio_minus_one`
- `variance_ratio_curve`

### Rank, sign, bootstrap, and panel procedures

- `wright_tests`
- `joint_wright_tests`
- `wright_critical_values`
- `joint_wright_critical_values`
- `variance_ratio_bootstrap`
- `automatic_vr_bootstrap`
- `subsample_variance_ratio`
- `panel_variance_ratio`

### Other martingale-difference tests

- `automatic_portmanteau`
- `average_exponential_test`
- `spectral_shape_test`
- `generalized_spectral_test`
- `dominguez_lobato_statistic`
- `dominguez_lobato_test`
- `chen_deo_test`

Bootstrap routines use Fortran's intrinsic random-number generator. Call
`seed_random(integer_seed)` before a bootstrap or simulation when repeatable
results are desired.

## Design choices

R lists are represented by named derived types. Array outputs are allocatable,
so callers are not required to choose maximum holding periods or bootstrap
sizes at compile time.

The original `VR.plot` numerical calculations are available through
`variance_ratio_curve`; graphics are intentionally not part of this numerical
library. The original `exrates.rda` file is retained for provenance but is not
read by the Fortran code.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for detailed mapping,
behavioral differences, and validation information.

## License

GPL version 2 only. See `LICENSE` and `NOTICE`. The original package metadata
and source are retained under `original/vrtest-1.2`.
