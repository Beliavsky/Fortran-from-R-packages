# sharpeRratio-fortran

A modern Fortran/FPM translation of the computational code in the R package
`sharpeRratio` 1.4.3.

The package estimates a signal-to-noise ratio without using sample moments.
It counts upper and lower records of cumulative returns across random
permutations, maps the normalized record imbalance through the package's
calibrated spline functions, and corrects for sample size and tail thickness.

## Implemented functionality

- upper and lower record counts;
- mean and confidence quantiles of the permuted record imbalance;
- the original `a`, `a_medium`, and `f` calibration splines;
- the Jelito-Pitera normality statistic used by the package;
- automatic Student-t tail-exponent fitting through the supplied `ghyp` port;
- fixed-tail and automatically fitted moment-free SNR estimation;
- filtering of NaN and infinite observations;
- deterministic seeded permutations;
- source-compatible and corrected confidence-quantile modes.

The public module is:

```fortran
use sharpe_rratio
```

All floating-point calculations use:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Basic example

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use sharpe_rratio

implicit none

real(dp) :: returns(100)
type(snr_result) :: fit
integer :: i

do i = 1, size(returns)
   returns(i) = 0.002_dp + 0.01_dp*sin(0.37_dp*real(i,dp))
end do

fit = estimate_snr(returns, num_perm=30, nu=5.0_dp, seed=17_int64)
if (.not. fit%ok) error stop trim(fit%message)

print *, fit%snr
print *, fit%ci_lower, fit%ci_upper
```

Omit `nu` to run the normality test and, when required, fit a skewed Student-t
model through `ghyp`:

```fortran
fit = estimate_snr(returns, num_perm=50, seed=42_int64)
```

## Main API

```fortran
integer = num_records_up(cumulative_values)
integer = num_records_down(cumulative_values)
records = compute_r0bar(x, num_perm, q1, q2, seed, source_compatible)
statistic = test_n(x)
value = a_full(r0)
value = f_full(x)
value = correction_b(r0,n)
value = theta_snr(r0,n,nu,nu_fixed)
nu = estimate_tail_exponent(x,...)
fit = estimate_snr(x,...)
```

`computeR0bar` and `estimateSNR` compatibility names are also provided.

## Source-compatible and corrected modes

By default, `estimate_snr` reproduces the upstream confidence-bound behavior:
the R argument named `quantiles` is ignored and the C++ defaults 0.025 and
0.975 are used with the original order-statistic indexing rule.

Set `source_compatible=.false.` to honor the supplied `quantiles` values and
use R type-7 interpolation. In corrected mode, omitted quantiles default to
0.05 and 0.95, matching the R function signature.

The random generator is not intended to reproduce `std::mt19937` draws from a
particular R session. Supplying `seed` gives reproducible Fortran results.

## Build

With FPM:

```text
fpm build
fpm test
fpm run sharpe_rratio_demo
fpm run --example known_tail_exponent
```

With GNU Fortran and Make:

```text
make test
make optimized
make demo
```

Windows users can run:

```text
scripts\validate.bat
```

## Scope

Rcpp, R lists, S3/R documentation machinery, and serialized R closures are not
required at runtime. The numerical calibration carried by the serialized
closures is compiled directly into the Fortran library. The upstream source
and data remain in `original/` for provenance.

See `API_MAP.md`, `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md`.

## License

GPL-3.0-only. See `LICENSE` and `NOTICE`.
