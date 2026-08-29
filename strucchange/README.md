# strucchange for modern Fortran

This package translates the computational core of R package **strucchange
1.6-0** to modern free-form Fortran and FPM. It provides structural-change
statistics, breakpoint dynamic programming, fluctuation processes, asymptotic
p-values and critical values, breakpoint confidence intervals, and monitoring
boundary calculations without requiring R.

The port preserves the upstream `GPL-2 | GPL-3` license choice and attribution.
See `NOTICE.md`, `UPSTREAM.md`, and `API_COVERAGE.md`.

## Dependencies

Place this directory at the repository root next to the shared packages:

```text
Fortran-from-R-packages/
  rfortran-core/
  rfortran-linalg/
  strucchange/
```

`fpm.toml` uses sibling path dependencies:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

No dependency source is vendored in this directory.

## Build and test

```sh
fpm build
fpm test
```

Run examples with:

```sh
fpm run --example breakpoints-example
fpm run --example cusum-example
```

The release hygiene script repeats the required build/test/style checks:

```sh
python tools/release_check.py
```

## Real kind

The package uses one real kind everywhere:

```fortran
use r_kinds, only : dp
```

The public `strucchange` module re-exports `dp`. Real variables use `real(dp)`
and real literals use `_dp` suffixes. There are no `double precision`,
`real*8`, `d0`/`D` exponents, or package-local competing real kinds.

## Small breakpoint example

```fortran
program demo
   use strucchange, only : dp, breakpoint_result, compute_breakpoints
   implicit none
   real(dp) :: x(20, 2), y(20)
   type(breakpoint_result) :: fit
   integer :: i

   do i = 1, 20
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i, dp)
      y(i) = 1.0_dp + 0.2_dp*real(i, dp)
      if (i > 10) y(i) = y(i) + 1.0_dp
   end do

   call compute_breakpoints(x, y, 4, 1, fit)
   if (fit%info == 0) print *, fit%breakpoints
end program demo
```

## Design

R-level formula/S3/plotting machinery is intentionally not recreated. APIs take
ordinary Fortran arrays and return arrays or small derived result types. Shared
statistics and linear algebra come from `rfortran-core` and `rfortran-linalg`.
See `API_COVERAGE.md` for the exact mapping and explicit exclusions.
