# Testing

## FPM

```text
fpm test
```

The test suite contains:

1. `test_connected_graph`
   - exact connected covariance recovery
   - convergence, support, symmetry, and zero row sums
2. `test_regular_heavytail`
   - Gaussian and Student-t simulated observations
   - residual and elapsed-time histories
3. `test_kcomp_heavytail`
   - Gaussian and Student-t two-component graphs
   - spectral component checks and optional objective history
4. `test_options_and_helpers`
   - QP initialization
   - scalar/vector option handling
   - invalid Student-t degrees of freedom
   - scalar weight helper
   - safe `maxiter` exhaustion

The k-component tests check the defining Laplacian spectrum rather than requiring
one exact edge support. This avoids the portability-sensitive assertion that was
identified in the earlier spectralGraphTopology translation.

## Strict GNU Fortran validation

```text
./run_gfortran_tests.sh
```

The script uses runtime checking, floating-point traps, warnings as errors, and a
separate optimized build. It runs all tests, the demo, and every example.

The Windows batch file performs the corresponding strict test build with
`gfortran`.
