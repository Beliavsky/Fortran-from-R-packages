# KFAS -- computational kernels in modern Fortran

This directory packages the computational engine of the R package **KFAS
1.6.0** for direct use from Fortran with the Fortran Package Manager (FPM).
KFAS implements Kalman filtering and smoothing for Gaussian and exponential-
family state-space models, including exact diffuse initialization and the
non-Gaussian approximation / importance-sampling machinery.

The upstream numerical engine was already written primarily in free-form
Fortran. This port preserves those algorithms, removes the R registration and
Rmath dependency, adds a Fortran-facing module, deterministic tests, examples,
and FPM packaging, and retains license/provenance material.

## Requirements

- A Fortran 2008-capable compiler (tested with GNU Fortran in the translation
  environment).
- FPM.

FPM obtains the pinned `fortran-lapack` dependency automatically. No system
BLAS/LAPACK installation or copied dependency source is required.

The manifest enables FPM's `implicit-external` compatibility switch because the
retained upstream numerical kernels call the classic BLAS/LAPACK and KFAS
external-procedure ABIs. A small compatibility layer delegates those calls to
`fortran-lapack`; it does not reimplement the numerical routines. The new
`kfas` public module itself uses explicit interfaces for its exposed wrappers.
All maintained numerical source uses the package-wide `dp = real64` kind from
`kfas_kinds`; the public `kfas` module re-exports `dp` for applications.

## Build and test

```text
fpm build
fpm test
```

Run the deterministic examples with:

```text
fpm run --example local-level
fpm run --example poisson-fixed
```

## Modern API example

```fortran
use kfas, only: dp, kfas_model, kfas_filter_result, kfas_gaussian_filter

type(kfas_model) :: model
type(kfas_filter_result) :: result
integer :: info

! Allocate and fill model%y, model%missing, model%z, model%h,
! model%tmat, model%rmat, model%q, model%a1, model%p1, model%p1inf.
call kfas_gaussian_filter(model, result, filter_signal = .true., info = info)
```

`kfas_gaussian_filter`, `kfas_gaussian_loglik`, and `kfas_gaussian_smooth`
automatically apply KFAS's LDL observation-equation transformation when a
multivariate `H` has nonzero off-diagonal terms. The module also exposes native
non-Gaussian Gaussian-approximation and Laplace-likelihood entry points for the
Poisson, binomial, gamma, and negative-binomial families; see
`example/poisson_fixed.f90` and `API_COVERAGE.md`.

Array conventions follow the upstream kernel conventions. For a model with
`n` time points, observation dimension `p`, state dimension `m`, and disturbance
dimension `r`:

- `y(n,p)` and `missing(n,p)`
- `z(p,m,nz)`
- `h(p,p,nh)`
- `tmat(m,m,nt)`
- `rmat(m,r,nr)`
- `q(r,r,nq)`
- `a1(m)`, `p1(m,m)`, `p1inf(m,m)`

For each system array the third dimension is 1 when constant and `n` when time
varying. Set the corresponding `time_varying` entry to 0 or 1. Set
`diffuse_rank` to the rank of `p1inf` when exact diffuse initialization is used.

See `example/local_level.f90`, `API_COVERAGE.md`, and `NOTICE.md`.

## License and citation

This adaptation follows the upstream `GPL (>= 2)` declaration and is licensed
under GPL-2.0-or-later. See `LICENSE`, `COPYING`, and `NOTICE.md`.

Please cite Jouni Helske (2017), "KFAS: Exponential Family State Space Models in
R," *Journal of Statistical Software*, 78(10), 1-39,
doi:10.18637/jss.v078.i10.
