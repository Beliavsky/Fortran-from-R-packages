# clarabel-fortran

Modern Fortran 2018 interface to the Clarabel conic interior-point solver, based on the computational API of the R package `clarabel` 0.11.2 and Clarabel.rs 0.11.1.

The attached R package is itself a language binding around the Rust solver. This port follows the same architecture:

- an idiomatic Fortran API for sparse problem data, cones, settings, results, one-shot solves, and persistent solver objects;
- an ISO C binding layer;
- a small Apache-2.0 Rust C-ABI bridge linked to the same vendored Clarabel.rs solver used by the upstream package.

It is therefore a full language binding, not a partial independent reimplementation of the interior-point method.

## Problem form

The solver minimizes

```text
0.5 * transpose(x) * P * x + transpose(q) * x
```

subject to

```text
A * x + s = b,  s in K.
```

`P` is supplied in upper-triangular CSC form. `K` can contain any ordered composition of:

- zero cones;
- nonnegative cones;
- second-order cones;
- exponential cones;
- power cones;
- generalized-power cones;
- positive-semidefinite triangle cones.

## Main Fortran API

```fortran
use clarabel

type(csc_matrix)              :: p, a
type(clarabel_cone)           :: cones(2)
type(clarabel_settings)       :: settings
type(clarabel_solution)       :: solution
type(clarabel_solver_type)    :: solver

p = csc_from_symmetric_upper(p_dense)
a = csc_from_dense(a_dense)
cones = [zero_cone(1), nonnegative_cone(4)]
settings = default_clarabel_settings()
settings%verbose = .false.

call clarabel_solve_problem(p, q, a, b, cones, solution, settings, code, message)
```

For sparse inputs, use `csc_from_arrays` or `csc_from_triplets`. Triplet construction sorts row indices, aggregates duplicates, and optionally drops small entries.

Persistent solvers support repeated solves and updates:

```fortran
settings%presolve_enable = .false.
settings%chordal_decomposition_enable = .false.
settings%input_sparse_dropzeros = .false.

call solver%initialize(p, q, a, b, cones, settings, code, message)
call solver%solve(solution, code, message)
call solver%update(q=q_new, b=b_new, code=code, message=message)
call solver%solve(solution, code, message)
```

Updates to `P` or `A` must preserve the original CSC sparsity pattern exactly.

## Building

### Fortran frontend

A normal FPM build now compiles and links without a prebuilt Rust library:

```sh
fpm build
```

The package contains `src/clarabel_dynamic_loader.c`, which loads the Clarabel
shared library only when a solver is created. This avoids the former
`cannot find -lclarabel_fortran_bridge` failure.

### Production Rust backend

A Rust toolchain, C compiler, BLAS, and LAPACK are required for actual Clarabel
solves. The dependency tree is bundled in `rust_bridge/vendor.tar.xz`, so the
Cargo build is offline and reproducible.

On Windows, the recommended first run is:

```bat
scripts\build_with_backend.bat run
```

This builds the production Rust DLL, copies any required GNU runtime DLLs beside
it, verifies the DLL through the Fortran interface, and then runs the demo.
After that one-time backend build, ordinary commands work:

```bat
fpm build
fpm test
fpm run
```

The wrapper can also invoke other FPM commands:

```bat
scripts\build_with_backend.bat build
scripts\build_with_backend.bat test
```

On GNU/Linux or macOS:

```sh
./scripts/build_backend.sh
fpm run
```

or:

```sh
./scripts/build_with_backend.sh build
./scripts/build_with_backend.sh run
./scripts/build_with_backend.sh test
```

The backend script copies the shared library to `rust_bridge/bin`. The runtime
loader checks that directory automatically when the current working directory
is the package root. For a program using this package as a dependency, set
`CLARABEL_FORTRAN_BRIDGE` to the absolute `.dll`, `.so`, or `.dylib` path, or
place the library on the normal operating-system search path.

`fpm build` and `fpm test` do not build or install the production Rust backend.
Without the shared backend, the default demo prints the exact setup command and
exits normally; library calls return `clarabel_backend_unavailable` (`-100`).
Passing frontend tests therefore does not imply that the production DLL exists.

If `rust_bridge\bin\clarabel_fortran_bridge.dll` exists but Windows still
reports that the module cannot be found, a dependent runtime DLL is missing.
Re-run `scripts\build_backend.bat`; the revised script copies common MinGW
runtime dependencies beside the bridge and verifies an actual QP solve.

## Frontend-only validation

A test-only C backend is included to validate the Fortran data model and ABI when Rust is unavailable:

```sh
./scripts/build_frontend.sh checked
./scripts/build_frontend.sh optimized
```

The mock solves only unconstrained and equality-constrained dense quadratic programs. It is never installed as the production solver and does not validate Clarabel's conic numerics.

## PSD vectorization

`psd_svec_upper` and `psd_smat_upper` implement Clarabel's scaled upper-triangular representation. Off-diagonal entries are multiplied by `sqrt(2)` in vector form.

## Directory layout

```text
src/             Fortran modules
include/         C ABI header
rust_bridge/     Rust bridge, lock file, and offline vendor archive
test/            frontend and ABI regression tests
integration/     production-backend numerical tests
example/         runnable examples
app/             demonstration program
scripts/         Unix and Windows build helpers
upstream/        original attached R package snapshot
docs/            API and porting documentation
```

## Scope differences from R

R S3/S4 objects, `Matrix` class dispatch, named lists, formula/data-frame handling, printing, and package documentation infrastructure are not reproduced. They are replaced by typed Fortran structures and explicit arrays. Cone ordering is always explicit in the Fortran `cones(:)` array, so no `strict_cone_order` flag is needed.

## License

The Fortran interface and Rust bridge are Apache-2.0. Clarabel.rs is Apache-2.0. The original R package snapshot and the vendored Rust dependencies retain their original notices and licenses.
