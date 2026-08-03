# osqp-fortran

Modern Fortran/FPM interface to the OSQP quadratic-programming solver, translated from the computational interface of the R package `osqp` 1.0.0.

OSQP solves convex quadratic programs of the form

```text
minimize    0.5*x'*P*x + q'*x
subject to  l <= A*x <= u
```

where `P` is positive semidefinite. Only the upper triangular part of `P` is stored and passed to OSQP.

## Features

- One-shot and persistent solver interfaces
- Dense and compressed-column sparse input
- Triplet-to-CSC conversion with sorting and duplicate aggregation
- Full OSQP 1.0 settings record
- Updates to `q`, `l`, `u`, and selected or complete `P`/`A` values
- Warm starts and cold starts
- Primal and dual solutions
- Primal- and dual-infeasibility certificates
- Solver status, residuals, duality gap, iteration counts, timing, and rho information
- Runtime backend loading, so the Fortran frontend builds without link flags
- Bundled OSQP 1.0.0, QDLDL, and AMD sources for offline backend builds

## Build

The Fortran frontend has no external link dependency:

```text
fpm build
fpm test
```

Tests that require the numerical backend report `SKIP` until the backend is built.

### Windows backend

From the package root:

```text
scripts\build_backend.bat
fpm run
```

or in one command:

```text
scripts\build_with_backend.bat run
```

The script uses CMake plus the same MinGW `gcc` toolchain as `gfortran`, builds the bundled OSQP source, creates `backend\bin\osqp_fortran_bridge.dll`, and verifies a real QP solve when FPM is available.

### Linux/macOS backend

```text
scripts/build_backend.sh
fpm run
```

The shared bridge is placed in `backend/lib`.

A bridge at another location can be selected with the environment variable `OSQP_FORTRAN_BRIDGE`.

## Basic example

```fortran
program example
   use osqp
   implicit none

   real(dp) :: p(2,2), q(2), a(3,2), l(3), u(3)
   type(osqp_solution) :: solution
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   p = reshape([4.0_dp, 1.0_dp, &
                1.0_dp, 2.0_dp], [2,2])
   q = [1.0_dp, 1.0_dp]
   a = reshape([1.0_dp, 1.0_dp, 0.0_dp, &
                1.0_dp, 0.0_dp, 1.0_dp], [3,2])
   l = [1.0_dp, 0.0_dp, 0.0_dp]
   u = [1.0_dp, 0.7_dp, 0.7_dp]

   settings%verbose = .false.
   call solve_osqp(q, solution, status, p=p, a=a, l=l, u=u, settings=settings)

   print *, trim(solution%status)
   print *, solution%x
end program example
```

See `example/` for persistent updates, warm starts, triplet sparse input, and infeasibility detection.

## Index conventions

The public Fortran sparse matrices use one-based row indices and one-based CSC column pointers. Matrix-update index arrays are also one-based. Conversion to OSQP's zero-based C convention occurs internally.

Matrices use ordinary Fortran `(row,column)` storage. Dense arrays are converted to CSC before entering the backend.

## Scope

The port covers the computational interface exposed by the R package: setup, solve, settings, updates, warm/cold starts, sparse coercion, solution information, and certificates. R S7/S3 objects, formula/data-frame handling, `Matrix` classes, R warnings, and deprecated `$` dispatch are omitted.

OSQP code generation and adjoint derivatives are disabled in the bundled backend because they are not exposed by the attached R package build. The backend uses OSQP's built-in direct QDLDL solver.

## License

Apache License 2.0. See `LICENSE`, `NOTICE.md`, and `THIRD_PARTY_LICENSES.md`.
