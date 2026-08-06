# Validation

The test suite contains five independent programs.

1. `block_top`: exact scalar fixture for the upstream block construction.
2. `conditional`: exact known covariance partition, regression slope,
   innovation scale, intercept, precision, and variance coefficients.
3. `static`: recovers an intercept near 1.2 and slope near 2.3 from a
   deterministic contaminated regression sample.
4. `dynamic`: exercises StAR, StDLM, and StVAR, including a two-column trend
   and all result dimensions.
5. `inference_errors`: checks numerical Hessian/Jacobian standard errors and
   invalid row/lag handling.

Both checked (`-fcheck=all`) and optimized (`-O3`) GNU Fortran 14.2 builds are
required to pass before packaging.
