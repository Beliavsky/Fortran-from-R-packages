# Testing

Four test programs cover:

1. Kernel formulas, string spectra, kernel operations, `sigest`, incomplete
   Cholesky, `ipop`, and probability coupling.
2. KPCA projection consistency, KCCA, kernel k-means, spectral clustering,
   ranking, and CSI.
3. LS-SVM regression/classification, SVM classification, Gaussian processes,
   RVM, kernel quantile regression, and online updates.
4. MMD, KHA, and KFA.

The scripts use GNU Fortran with standard conformance, warnings as errors,
runtime bounds checking, and floating-point traps. An optimized script repeats
all tests, examples, and the demo with `-O3 -Werror`.
