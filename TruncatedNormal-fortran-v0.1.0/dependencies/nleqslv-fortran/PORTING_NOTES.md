# Porting notes

## Numerical core

`nleqslv` already contains its principal nonlinear-equation algorithms in
Fortran. Those routines were retained algorithmically and converted from
fixed-form `.f` to standard free-form `.f90`. The R/C callback, allocation,
and object-construction layer was replaced by the modern module
`nleqslv_fortran`.

The original routines still intentionally retain their mature procedural
layout and BLAS/LAPACK calls. In v0.1.1 the converted kernels are collected
inside the module `nleqslv_legacy_kernels`, so calls among solver routines
have explicit interfaces. BLAS/LAPACK calls use the explicit interface module
`nleqslv_blas_lapack`. This is required for strict builds with
`-Werror=implicit-interface`; v0.1.0 incorrectly left those calls as legacy
implicit interfaces. The new public API uses `implicit none`, `real64`,
allocatable arrays, derived types, abstract procedure interfaces, and named
constants.

## Dependencies

The numerical kernels use BLAS and LAPACK (`DGEQRF`, `DORMQR`, `DORGQR`,
`DTRCON`, level-1/2/3 BLAS, etc.). FPM links system `lapack` and `blas` rather
than copying those libraries.

## Callbacks and reentrancy

The upstream solver callback ABI does not carry a user context pointer. The
modern wrapper therefore stores the currently active Fortran function and
Jacobian procedure pointers in module state while a solve is running. As in
the R package, recursive use is rejected. Concurrent calls from multiple
threads are not supported by this wrapper.

Using an internal procedure as a callback may cause some GNU toolchains to
request an executable stack for the caller's trampoline. Module or external
procedures avoid that toolchain behavior; the library objects themselves do
not require such a trampoline.

## Non-finite values

The initial point and initial function value are required to be finite. During
ordinary globalization, a non-finite function value is replaced by a very
large finite value so that the original solver backtracks. During finite-
difference Jacobian evaluation, a non-finite callback result is treated as a
fatal callback error, matching upstream intent.

## `searchZeros`

`search_zeros` follows the computational behavior of the R helper: only
termination-code-1 solutions are retained, solutions are de-duplicated after
decimal rounding, and the unrounded representatives are returned in
lexicographic order. R's `try`/condition bookkeeping has no direct Fortran
exception analogue, so fatal callback errors terminate rather than being
collected as `idxfatal` entries.

## `testnslv`

`test_nleqslv` evaluates requested combinations of Newton/Broyden and global
strategies and returns numerical counters, criterion values, termination
codes, and optional CPU timings. R dataframe construction and printing are
not reproduced.

## Iteration output

The R package formats elaborate iteration tables in C. The Fortran port
retains the underlying solver trace hooks but emits compact text only. This is
presentation behavior and does not alter the numerical algorithm.

## Termination codes

The current R-facing meanings are retained:

1. function criterion near zero
2. x-values within `xtol`
3. no better point found / stalled
4. iteration limit exceeded
5. Jacobian too ill-conditioned
6. Jacobian singular
7. Jacobian unusable
-10. supplied Jacobian probably incorrect

## Validation

The test suite covers:

- Dennis-Schnabel example with numerical and analytic Jacobians
- every Newton/Broyden globalization strategy
- banded finite-difference Jacobians
- singular-Jacobian handling
- multiple-root search
- the `testnslv`-style method/global sweep

Both runtime-checked and optimized builds are exercised during release
validation.
