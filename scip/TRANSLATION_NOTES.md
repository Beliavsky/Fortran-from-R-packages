# Translation notes

## Scope

`scip` 1.10.0-3 is an R interface to the SCIP Optimization Suite and vendors
SCIP 10.0.2 plus SoPlex 8.0.2.  Rewriting the whole SCIP framework in Fortran
would be a separate multi-year solver project, not a translation of the
package-owned interface layer.  This release therefore translates the full
R-exposed computational API to modern Fortran while compiling and linking the
same bundled SCIP/SoPlex source as the backend.

No alternative optimizer is substituted.  Linear, mixed-integer, nonlinear
quadratic, SOS, indicator, presolve, branching, heuristics, and solution-pool
behavior still comes from the exact solver sources bundled with the upload.

## R to Fortran mapping

- R external pointers -> `type(scip_model_t)` containing an opaque C pointer.
- `.Call` functions -> explicit `bind(C)` interfaces in `scip_c_api.f90`.
- R dense / `dgCMatrix` / `simple_triplet_matrix` conversion -> native dense
  arrays plus `type(scip_csc_matrix)` and `make_csc_matrix()`.
- `scip_solve()` -> generic Fortran `scip_solve()` for dense or CSC matrices.
- `scip_model()` -> constructor returning `scip_model_t`.
- `scip_add_*`, `scip_set_*`, `scip_get_*` -> both package-style module
  procedures and type-bound methods.
- `scip_control()` -> `type(scip_control)` with matching common fields and
  `control%add_param(name,value)` for additional native SCIP parameters.

## Memory ownership

A Fortran model owns an external SCIP allocation and must be released with
`call model%free()` / `call scip_model_free(model)`.  The type intentionally
has no automatic finalizer: intrinsic assignment of a derived type containing
an opaque pointer can otherwise free a temporary constructor result and leave
the assigned object dangling.  Treat `scip_model_t` as a non-copyable handle.
The one-shot `scip_solve()` routines release their internal model themselves.

## Backend build

The uploaded source is the CRAN-style R package source, where some upstream
SCIP/SoPlex source was patched to route output through R's printing functions
and several CMake-only directories were stripped.  `scripts/build_vendor.*`
handles both facts without requiring R:

1. temporary CMake stubs are created for stripped non-library directories;
2. SCIP and SoPlex are built as static libraries;
3. `standalone_streams.cpp` provides standard-output replacements for the
   small set of R printing/error symbols referenced by the patched sources;
4. the solver libraries and `scip_fortran_shim.c` are merged into
   `libscipfortran_backend.a`.

The original source tree itself remains under `original/`.

## Intentional interface differences

- R S3 printing and Matrix/slam class dispatch are omitted.
- Variable and constraint names are ordinary Fortran strings.
- Model memory is explicitly freed rather than garbage-collected.
- Fortran one-shot bounds/types are optional vectors.  Scalar recycling from R
  is represented by simply omitting defaults or constructing constant arrays.
- Native SCIP parameters are strongly typed through overloaded Fortran
  procedures instead of R's dynamically typed lists.

## Numerical behavior

The translation does not reimplement SCIP numerics.  All solver algorithms,
presolve, branching, cuts, heuristics, LP solves, and nonlinear handling are
provided by the bundled SCIP 10.0.2 and SoPlex 8.0.2 sources.
