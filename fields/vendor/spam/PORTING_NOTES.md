# Porting notes

## Translation strategy

All Fortran numerical sources shipped in `spam/src` were converted from fixed form to standard free form and retained under `src/native/`. Comments and routine names are preserved. The conversion changes statement layout/continuations only, except for one fixed-form blank-insensitive numeric literal (`2 000 000 000`) which necessarily becomes `2000000000` in free form. `tools/convert_fixed.py` is included for auditability.

Modern modules then provide typed APIs over the native routines and R-level computational formulas. The R/C++ registration layer is not required. The Rcpp nearest-distance implementation is not separately ported because the package already ships the equivalent Fortran distance kernel and the modern wrapper implements the same four metrics directly without fixed output-buffer guessing.

## Sparse representation

`csr_matrix` stores 1-based CSR arrays, matching spam's `entries`, `colindices`, and `rowpointers` conventions. Triplet conversion sorts columns within rows, consolidates duplicate `(i,j)` entries by summation, and can deliberately retain stored numerical zeros when `eps < 0`. This is important for `nearest.dist`, because spam intentionally preserves zero distances such as diagonal entries so subsequent covariance transforms can turn them into nonzero variances.

The public implementation uses default Fortran integers. It therefore corresponds to the ordinary `spam` package; the separate `spam64` extension and dotCall64 64-bit storage-selection machinery are outside this port.

## Cholesky

Ng-Peyton/PCx factorization is not replaced by dense LAPACK. `spam_chol_factor` calls the translated `cholstepwise` routine and dynamically retries if the factor/index work arrays are too small. `spam_chol_update` calls the upstream `updatefactor`; solve routines call the original block triangular kernels. `spam_logdet` returns the log determinant of the original SPD matrix, i.e. twice the log determinant of the Cholesky factor. `chol_factor_logdet` corresponds to the upstream determinant of a `spam.chol.NgPeyton` factor itself.

## ARPACK

The complete symmetric and nonsymmetric ARPACK sources bundled by spam are retained. `spam_eigen_symmetric` and `spam_eigen_general` wrap the package's own `ds_eigen_f`/`dn_eigen_f` drivers and expose convergence counts/iterations through `eigen_result`.

## Random generation

Covariance-form Normal simulation uses a dense LAPACK Cholesky after converting the supplied covariance matrix to dense form. Precision-form and canonical simulation retain sparse Ng-Peyton factorization/triangular solves. This mirrors the underlying algorithms but means covariance-form sampling is not the scalable sparse path; use precision-form simulation when sparse scalability matters.

The constrained samplers use the standard Rue-Held projection identity. Conditional simulation uses the exact Gaussian conditional mean/covariance identity. `rgrf` supports the Cholesky method, which is also the only method implemented in the attached upstream `rgrf.R`.

## MLE

The Gaussian likelihood is the upstream expression

`n*log(2*pi) + log(det(Sigma)) + r' Sigma^{-1} r`.

The Fortran MLE exposes both the built-in spam covariance families by name and custom typed covariance callbacks (`mle_spam_custom`). It uses a bounded finite-difference BFGS implementation because the supplied `r_mod` optimization callback interface is declared `pure`, while sparse Cholesky factorization is necessarily impure. This is a genuinely missing compatibility layer rather than a reimplementation of probability/statistical helpers already provided by `r_mod`.

## R interfaces intentionally not reproduced

S4/S3 object construction, formula/model-frame behavior, plotting/grid/display functions, Matrix-class coercions, print/summary methods, dynamic R option management (`powerboost`), package-version helpers, and foreign MatrixMarket/Harwell-Boeing IO wrappers are interface/presentation concerns and are omitted. Generic `apply.spam` is represented by typed Fortran callbacks; arbitrary dynamically shaped R return values are not emulated.
