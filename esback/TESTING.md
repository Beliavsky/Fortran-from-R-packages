# Testing

Validated with GNU Fortran 14.2 in two clean configurations:

- strict: warnings as errors, bounds/runtime checks, and invalid/zero/overflow traps
- optimized: `-O3`, warnings as errors

Permanent tests:

- `test_er_cc`: fixed independent Python references for both exceedance-residual variants and all four conditional-calibration p-values
- `test_esreg_core`: specification functions and a fixed Fissler-Ziegel loss reference
- `test_esr`: strict, auxiliary, and intercept ESR fits; coefficient references from an independent SciPy optimization; covariance and p-value checks
- `test_nuisance`: conditional location-scale fitting, quantile regression, NID sparsity, empirical conditional CDF, and truncated-normal variance

Run:

```sh
./tools/test_gfortran.sh strict
./tools/test_gfortran.sh optimized
```
