# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational-core translation of `bigstatsr` 1.6.2.
- Added portable stream-backed real64 FBMs and 0..255 coded FBMs.
- Added blocked products, crossproducts, correlation, scaling, counts and transpose.
- Added AUC/bootstrap AUC and partial correlation.
- Added univariate linear and logistic regression.
- Added Gaussian and logistic elastic-net/lasso paths.
- Added primal/dual partial SVD and matrix-free ARPACK SVD via vendored RSpectra Fortran.
- Added `get_beta`, `block_size`, and grouped sparse-regression summaries.
- Added strict numerical tests and FPM metadata.
