# mice - modern Fortran computational translation

This directory translates reusable numerical algorithms from R package
**mice 3.19.0** to modern free-form Fortran with FPM. It is intended to live as
a top-level directory beside the shared packages in
`Fortran-from-R-packages`.

The public module is `mice`; it re-exports `dp` and the principal numerical
routines. Missing numeric cells are represented by IEEE NaNs where a complete
matrix is supplied. Lower-level imputers also accept explicit logical
observed/where masks.

The current numerical surface includes the FCS driver, normal and PMM families,
`midastouch`, `mpmm`, binary/multinomial/proportional-odds/LDA categorical
imputation, NARFCS offset kernels, the heterogeneous `2l.norm` Gibbs sampler,
`2lonly.*`, pooling, missingness diagnostics, amputation, and native matching
and Legendre helpers. See `API_COVERAGE.md` for exact scope and remaining
external-model gaps.

## Build

From this directory in the target repository:

```text
fpm build
fpm test
```

The FPM manifest uses sibling dependencies only:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

No system BLAS/LAPACK installation or `-lblas`/`-llapack` link flag is used by
this package. `rfortran-linalg` supplies the shared pure-Fortran LAPACK-backed
linear algebra used here.

## FCS example

```fortran
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
use iso_fortran_env, only : int64
use mice, only : dp, mice_fcs_result, mice_fcs_impute, &
                 mice_method_pmm, mice_method_norm

real(dp) :: data(6, 2), nan
integer :: methods(2), predictors(2, 2), info
type(mice_fcs_result) :: result

nan = ieee_value(0.0_dp, ieee_quiet_nan)
data(:, 1) = [1.0_dp, 2.0_dp, nan, 4.0_dp, 5.0_dp, 6.0_dp]
data(:, 2) = [2.0_dp, nan, 6.0_dp, 8.0_dp, 10.0_dp, 12.0_dp]
methods = [mice_method_pmm, mice_method_norm]
predictors = reshape([0, 1, 1, 0], [2, 2])
call mice_fcs_impute(data, methods, predictors, 5, 2, 1234_int64, result, info)
```

See `example/mice_example.f90` and `example/mice_parity_example.f90` for
complete programs.

## Licensing and provenance

The upstream package is GPL (>= 2). This translation is GPL-2.0-or-later.
See `LICENSE`, `NOTICE.md`, `PROVENANCE.md`, and the retained upstream metadata
and computational reference files under `upstream/`.
