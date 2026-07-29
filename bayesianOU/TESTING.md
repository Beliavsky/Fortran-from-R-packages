# Testing

## Commands

```text
./build_gfortran.sh strict
./build_gfortran.sh release
```

Strict flags include:

```text
-std=f2018 -Wall -Wextra -Werror -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

Release uses `-O3` with warnings treated as errors.

## Test programs

- `test_utils`: training-window standardization, PCA factor construction,
  orthogonalization, and level specifications.
- `test_geometry`: Euclidean HMC on an anisotropic Gaussian plus SoftAbs mass
  construction.
- `test_reference`: independently calculated fixed Student-t OU log-likelihood
  and pointwise values.
- `test_model`: single-level simulation, fit, posterior summaries, PSIS-style
  LOO, and structural validation.
- `test_nested`: full three-level simulation and fit, Level-2 parameters, latent
  mean trajectory, and OOS output.
- `test_mi`: multiple-imputation fitting and Rubin pooling.
- `test_psis`: deterministic PSIS-style LOO smoke/reference checks.

The build script also compiles and runs:

- `demo_bayesianou`
- `nested_three_level`
- `geometry_hmc`

## Results

Both strict and optimized builds pass all tests and runnable targets with GNU
Fortran 14.2 and the system BLAS/LAPACK libraries.

The original R package relies on Stan and optional R packages that are not
installed or needed for the Fortran tests. The original package tree is retained
under `original/` for source-level comparison.
