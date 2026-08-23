# RSpectra-fortran

A modern Fortran/FPM translation of the computational interface of CRAN `RSpectra` 0.16-2.

`RSpectra` is an R interface to the C++ Spectra library for computing a small number of eigenvalues/eigenvectors or singular values/vectors of large matrices. This port keeps the same computational model but uses the mature Fortran ARPACK reverse-communication solvers as the iterative backend and LAPACK/BLAS for factorizations and dense full-spectrum fallbacks.

## Implemented API

The top-level module is:

```fortran
use rspectra
```

It exports:

- `eigs_sym(A, k, which, sigma, opts)` for real symmetric problems
- `eigs(A, k, which, sigma, opts)` for real general problems with complex-valued results
- `svds(A, k, nu, nv, opts, ...)` for truncated SVD
- dense matrix, CSR sparse matrix, and user-defined operator interfaces
- `eigs_opts` and `svds_opts`
- result types carrying values/vectors, `nconv`, `niter`, `nops`, and `info`

### Eigenvalue selection rules

Symmetric problems support the RSpectra rules:

- `LM`: largest magnitude
- `SM`: smallest magnitude
- `LA`: largest algebraic
- `SA`: smallest algebraic
- `BE`: both ends, ordered largest, smallest, second largest, second smallest, ...

General problems support:

- `LM`, `SM`
- `LR`, `SR`
- `LI`, `SI`, where imaginary part is compared in magnitude, matching Spectra

### Shift-and-invert

For symmetric dense matrices, real shift-and-invert is implemented with a LAPACK LU factorization of `A - sigma I` and ARPACK mode 3.

For CSR symmetric matrices, the shifted path currently converts the matrix to dense form before factorization. Ordinary unshifted CSR problems remain sparse throughout.

For general explicit matrices, shifted selection is computed by a full LAPACK eigen-decomposition followed by exact sorting according to `1 / (lambda - sigma)`. This supports both real and complex shifts and gives the intended RSpectra selection semantics, but it is an O(n^3) fallback rather than a large-scale iterative shifted solve.

Shift mode is intentionally not offered for an arbitrary user-defined operator because RSpectra itself requires a solver for `(A - sigma I)x=b`, not merely an `A*x` callback.

### SVD

`svds` follows the same normal-equation strategy as RSpectra:

- tall matrices use `B'B`
- wide matrices use `BB'`
- the other singular-vector side is recovered by a matrix operation and division by the singular value
- optional centering/scaling uses `B = (A - 1 c') diag(1/s)`

Dense `center=.true.` computes column means. Dense `scale=.true.` computes column Euclidean norms after centering, matching RSpectra. For abstract operators, provide explicit center/scale vectors.

## User-defined matrix operators

Instead of an R function object, extend the abstract type `linear_operator` and implement `prod` and `tprod`:

```fortran
type, extends(linear_operator) :: my_operator
    ! user data
contains
    procedure :: prod => my_prod
    procedure :: tprod => my_tprod
end type my_operator
```

Then call `eigs`, `eigs_sym`, or `svds` directly on the operator instance.

## CSR format

`make_csr_operator()` accepts conventional Fortran 1-based CSR arrays:

- `row_ptr(1) = 1`
- row `i` occupies entries `row_ptr(i):row_ptr(i+1)-1`
- `col_ind` is 1-based

## Build requirements

The FPM manifest links these system libraries:

```toml
link = ["arpack", "lapack", "blas"]
```

On a typical Linux system this corresponds to ARPACK-NG plus LAPACK/BLAS development packages.

Build and test with FPM:

```text
fpm test
```

or compile directly with GNU Fortran and link `-larpack -llapack -lblas`.

## Validation

The supplied tests cover:

- symmetric `LM`, `SM`, and `BE`
- general real and complex eigenpairs
- real symmetric shift-invert
- general complex-shift selection
- CSR sparse matrices
- user-defined operator interface
- tall and wide truncated SVD
- centering/scaling
- full-spectrum LAPACK fallback
- eigenvector and singular-vector residuals

The project is tested with strict GNU Fortran 2018 checks.

## Scope differences from R

R S3 dispatch and Matrix-package classes are not translated. The equivalent numerical interfaces are dense arrays, CSR arrays, and Fortran operator types. Automatic R class coercion, warnings, and printed list objects are presentation/interface behavior rather than numerical algorithms and are omitted.

See `PORTING_NOTES.md` and `docs/TRANSLATION_COVERAGE.md` for details.
