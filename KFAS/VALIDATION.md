# Validation record

The original translation validation used GNU Fortran with system BLAS/LAPACK.
The package has since been built and tested with FPM on Windows using the pinned
`fortran-lapack` dependency, without system BLAS/LAPACK linker flags.

The following equivalent validation was completed directly with GNU Fortran:

- every `src/*.f90` unit compiled;
- the complete library linked against LAPACK and BLAS;
- `test/test_kfas.f90` linked and ran successfully;
- `example/local_level.f90` linked and ran successfully;
- the Gaussian reference log-likelihood was `-7.275349449509625`;
- the exact-diffuse reference test passed;
- Gaussian smoothing matched an independently computed RTS reference for state
  means and variances;
- correlated multivariate `H` matched independently computed multivariate
  Gaussian log-likelihoods both with complete data and with a partially missing
  observation after the native LDL transformation;
- a fixed-state Poisson model matched its exact analytic Poisson likelihood
  through the non-Gaussian approximation/log-likelihood wrappers;
- the native normal, Poisson, binomial, gamma, and negative-binomial density
  bridge tests passed;
- `fpm build`, `fpm test`, and `fpm run --example --all` completed successfully
  after replacing system BLAS/LAPACK linkage with `fortran-lapack`;
- the same checks completed after standardizing maintained numerical source on
  the package-wide `dp = real64` kind;
- no copied dependency sources, duplicate Fortran files, semicolon-containing
  Fortran lines, over-132-column Fortran lines, fixed-form files, build products,
  caches, or nested ZIP files were present.

`tools/release_check.py` performs the requested literal `fpm build`, `fpm test`,
hygiene checks, and `fpm clean --all` on any machine with FPM installed.
