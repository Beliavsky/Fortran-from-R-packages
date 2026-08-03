# Porting notes

## Why the solver core remains Rust

The R package does not implement Clarabel's numerical algorithm in R. It validates and converts R objects and invokes Clarabel.rs. Replacing that architecture with a new Fortran interior-point solver would not be a translation of the attached package and would be much harder to verify. This port instead replaces the R-specific layer with a Fortran-specific layer and calls the same Apache-licensed solver core through a narrow C ABI.

## Matrix orientation and indexing

Fortran dense arrays are passed in conventional `(row,column)` form. CSC row indices and column pointers stored in `csc_matrix` are zero-based because that is the native Clarabel representation. `csc_from_arrays` can convert either zero- or one-based integer input. `csc_from_triplets` defaults to one-based Fortran indices.

`P` must be square and upper triangular in CSC form. `csc_from_symmetric_upper` checks symmetry and extracts the upper triangle. Clarabel ignores no implied lower-triangular entries: the upper triangle is the authoritative quadratic matrix representation.

## Cone order

The R package can sanitize a named cone list into a conventional order or pass repeated cones in explicit order. Fortran has no named-list ambiguity: `cones(1), cones(2), ...` correspond directly to consecutive row blocks of `A` and `b`.

## Settings defaults

The Fortran defaults follow `clarabel_control()` in the attached R package, including QDLDL selection and chordal decomposition disabled by default. `time_limit` uses the largest finite `real(dp)` value as the portable Fortran equivalent of R's `Inf`.

The low-level Rust solver has a few different native defaults. Every setting field is transmitted explicitly, so the Fortran behavior follows the R wrapper rather than relying on Rust defaults.

## Persistent updates

The Rust solver updates numeric values while preserving its internal sparse structure. The Fortran object stores the original CSC index arrays and rejects a `P` or `A` update if dimensions, nonzero count, column pointers, or row indices differ. Set all three of the following to false before construction when updates must be guaranteed:

```fortran
settings%presolve_enable = .false.
settings%chordal_decomposition_enable = .false.
settings%input_sparse_dropzeros = .false.
```

## Error handling

The Rust bridge catches Rust panics and converts errors to integer return codes plus an error string. The high-level Fortran API reports these through optional `code` and allocatable `message` arguments. Array and cone dimensions are validated before crossing the ABI.

## Validation limits in the translation environment

GNU Fortran and GCC were available, but Cargo/Rust was not. Consequently:

- all Fortran modules, C ABI declarations, examples, and production-backend integration sources were compiled;
- seven regression programs were executed in checked and optimized modes against the test-only ABI backend;
- the actual Clarabel.rs bridge could not be compiled or numerically executed in this environment.

The archive includes the offline vendor tree, a locked Cargo manifest, and production integration tests adapted from the upstream R tests so that a machine with Rust can perform the full validation using `scripts/build_with_backend.sh`.
