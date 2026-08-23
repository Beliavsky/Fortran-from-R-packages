# Porting notes

## Numerical backend

Upstream RSpectra wraps the C++ Spectra template library. A direct Fortran call into Spectra would require retaining a C++ template/shim layer and would not meaningfully translate the computational interface to Fortran. This port instead maps the same large-scale solver semantics onto ARPACK's native Fortran reverse-communication interface.

This is not a line-by-line translation of Spectra internals. It is a translation of RSpectra's computational contract to an established Fortran eigensolver backend.

## Symmetric eigenproblems

Regular `eigs_sym` uses `DSAUPD/DSEUPD`. Real shift-and-invert uses ARPACK mode 3 with a LAPACK factorization of `A - sigma I`.

The RSpectra `BE` ordering is preserved as:

```text
largest, smallest, second largest, second smallest, ...
```

## General eigenproblems

Regular `eigs` uses `DNAUPD/DNEUPD` and converts ARPACK's real storage of conjugate eigenvectors into Fortran complex vectors.

`LI` and `SI` compare `abs(Im(lambda))`, matching Spectra's `LARGEST_IMAG` and `SMALLEST_IMAG` selection rules.

### Shifted general matrices

ARPACK's real-nonsymmetric complex-shift modes use real or imaginary parts of the complex resolvent and do not reproduce Spectra's complex-shift ordering/result extraction directly. During validation this produced transformed rather than original eigenvalues for a simple rotation block.

For correctness, explicit general shifted matrices therefore use LAPACK `DGEEV` to compute the spectrum and then apply the exact RSpectra criterion to:

```text
1 / (lambda - sigma)
```

This makes real and complex shifts correct but sacrifices the large-scale advantage in that particular path. CSR matrices are densified for general shift mode.

## SVD

The implementation follows upstream `svds.cpp` conceptually. It solves a symmetric eigenproblem on `B'B` or `BB'`, where:

```text
B = (A - 1 c') diag(1/s)
```

and computes the opposite side with one additional matrix operation per requested vector.

## Sparse storage

R's `dgCMatrix`, `dgRMatrix`, `dsCMatrix`, and `dsRMatrix` S3 methods are replaced by a native CSR operator. Unshifted iterative operations stay sparse. Shifted CSR currently densifies because this v0.1 port does not bundle a sparse direct solver.

## Full-spectrum requests

RSpectra delegates requests for all or almost all eigenvalues/singular values to R's dense `eigen()`/`svd()`. The Fortran port similarly delegates to LAPACK `DSYEV`, `DGEEV`, and `DGESVD`.

## Initial vectors and diagnostics

`eigs_opts%initvec` maps to ARPACK's user-supplied residual/start vector. `niter` and `nops` are read from ARPACK `IPARAM` counters.

## Numerical differences

Eigenvectors and singular vectors are only unique up to sign or complex phase. Tests therefore emphasize residuals and eigenvalue/singular-value agreement rather than raw vector signs.

ARPACK and Spectra use different restart implementations, so iteration counts and operation counts need not be identical even when converged values agree.
