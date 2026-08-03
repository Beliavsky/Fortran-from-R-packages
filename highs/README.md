# highs-fortran

A modern Fortran/FPM interface to HiGHS 1.14.0, translated from the computational API of the R package `highs` 1.14.0-2.

The R package is a wrapper around the external HiGHS solver. This project follows the same architecture:

- idiomatic Fortran model, sparse-matrix, solution, basis, and solver types;
- a small C runtime loader compiled by FPM;
- a C++ bridge to the bundled HiGHS source;
- no direct native-library linker dependency during `fpm build`.

## Build the Fortran frontend

```text
fpm build
fpm test
```

These commands do not compile HiGHS. Solver-dependent tests print `SKIP` when the runtime backend is absent; sparse/model tests still run.

## Build the numerical backend

The source archive already contains HiGHS 1.14.0, so no network download is required.

### Windows with MinGW/GNU Fortran

```bat
scripts\build_backend.bat
fpm run
```

Or:

```bat
scripts\build_with_backend.bat run
```

The script requires CMake, GCC, and G++. It builds static HiGHS, creates `backend\bin\highs_fortran_bridge.dll`, copies common MinGW runtime DLLs, and verifies a real LP solve when FPM is available.

### Unix-like systems

```sh
scripts/build_backend.sh
fpm run
```

The bridge is placed in `backend/lib`. Set `HIGHS_FORTRAN_BRIDGE` to an explicit shared-library path when installing it elsewhere.

## High-level solve

```fortran
use highs
implicit none

type(highs_solution) :: solution
integer(highs_int) :: status
real(dp) :: a(1,2)

a = reshape([2.0_dp, 1.0_dp], [1,2])
call highs_solve([3.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], &
   solution, status, a=a, lhs=[-highs_default_infinity], rhs=[2.0_dp], &
   vartype=[highs_var_integer, highs_var_integer], maximum=.true.)
```

`highs_solve` supports:

- linear objectives;
- quadratic objectives using the HiGHS convention `0.5*x'*Q*x + c'*x + offset`;
- ranged linear constraints `lhs <= A*x <= rhs`;
- continuous, integer, semi-continuous, semi-integer, and implicit-integer variables;
- warm-start column values and solver controls.

## Persistent solver API

Use `highs_new_solver`, `highs_pass_model`, `highs_run`, and `highs_get_solution` when repeated solves or hot starts are needed. The persistent interface also supports:

- changing objective coefficients, variable bounds, row bounds, matrix coefficients, integrality, sense, and offset;
- setting typed HiGHS options;
- basis get/set/clear;
- primal and dual rays;
- presolve;
- model read/write;
- primal/dual start values.

Fortran-facing indices are 1-based. Sparse `start` and `index` arrays stored in `highs_sparse_matrix` use the 0-based HiGHS representation.

## Matrix orientation

Dense matrices use ordinary Fortran shape `(number_of_rows, number_of_columns)`. `highs_csc_from_dense`, `highs_csr_from_dense`, and `highs_csc_from_triplets` create canonical sparse storage and combine duplicate triplets.

## Backend discovery

The loader checks:

1. an explicit path passed to `highs_load_backend`;
2. `HIGHS_FORTRAN_BRIDGE`;
3. the package-root `backend/bin` or `backend/lib` locations;
4. the platform library search path.

A missing backend does not prevent compilation. `demo_highs` prints the setup command and exits normally.

## Validation

The delivered tests cover sparse conversion, model validation, backend/version loading, LP, QP, MIP, persistent model edits, basis transfer, and model-file round trips. The numerical tests were run against the bundled HiGHS 1.14.0 backend with GNU Fortran and GNU C/C++.

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for detailed coverage and design differences.
