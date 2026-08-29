# Validation

Validation environment:

- GNU Fortran 14.2.0
- Fortran 2018
- runtime checking enabled (`-fcheck=all`)
- implicit external interfaces treated as errors (`-Werror=implicit-interface`)
- BLAS/LAPACK linked

Passing retained programs:

1. `test_core`
2. `test_dependence`
3. `test_pairwise`
4. `test_copula_kriging`
5. `test_sim_fit`
6. `test_extended`
7. `test_parity_targets`
8. `test_remaining_gaps`
9. `test_conditional_full`
10. `test_design_inference`
11. `example/demo.f90`

Version 0.2.0 adds regression/smoke checks for:

- set-partition enumeration (`Bell(4) = 15`);
- conditional max-linear simulations exactly honoring conditioning values in a deterministic design case;
- Van der Corput unit-vector construction;
- 2D point and Cartesian-grid turning-band Gaussian simulation;
- Schlather, geometric-Gaussian, and extremal-t TBM max-stable paths;
- generic sandwich covariance/SE construction;
- Brown-Resnick stationary sandwich standard errors;
- equality between the sum of full-margin pair contributions and the aggregate Brown-Resnick log likelihood;
- latent DIC identities for an identical chain.

The supplied `r_mod.F90` causes a linker warning about an executable stack in this toolchain; it predates the new parity modules and does not cause a test failure. FPM was not installed in the validation container, so the same source dependency graph was compiled directly with `gfortran`; `fpm.toml` was retained for normal FPM builds and parsed as TOML during packaging.


## Version 0.3.0 gap tests

`test_remaining_gaps.f90` exercises the new latent-GEV MCMC driver, the FFT/circulant Gaussian-grid sampler (including an empirical marginal-variance regression check), and Schlather conditional simulation for a supplied hitting partition. All eight test programs pass under GNU Fortran 14.2.0 with `-std=f2018 -fcheck=all -Werror=implicit-interface`, linked with BLAS/LAPACK.

## Version 0.4.0 parity-closure tests

`test_conditional_full.f90` checks the translated multivariate Normal/Student probability kernels, normalized Schlather/extremal-t/Brown-Resnick partition weights, Gibbs partition updates, automatic conditional simulation for all three model families, starting hitting partitions, and the exact Brown-Resnick hitting-scenario output. Conditioning values are checked for exact preservation.

`test_design_inference.f90` verifies spatial-plus-temporal GEV design expansion by round-tripping known unit-Frechet values, exercises active-parameter sandwich inference for Smith, Schlather, Schlather-independence, Brown-Resnick, geometric-Gaussian, extremal-t, and spatial GEV, and checks generalized-Cauchy `smooth2` inference.

A clean direct `gfortran` build of all modules passed all ten test programs and `example/demo.f90` under GNU Fortran 14.2.0 with `-std=f2018 -fcheck=all -Werror=implicit-interface`, linked with BLAS/LAPACK. FPM is not installed in this environment; `fpm.toml` is retained and parsed as TOML during packaging.
