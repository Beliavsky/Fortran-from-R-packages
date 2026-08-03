# Testing

The test suite contains five programs.

1. `test_cla_reference`
   - Reproduces the documented TLT/VTI/GLD three-asset example.
   - Checks all turning-point weights, lambdas, means, and standard deviations
     against fixed independent values.
   - Checks budget, box bounds, and decreasing lambdas.
2. `test_cla_queries`
   - Checks mean-to-risk interpolation and risk-to-mean inversion.
3. `test_cla_bounds`
   - Exercises nonuniform lower/upper bounds and infeasibility detection.
4. `test_cla_garch`
   - Exercises normal and standardized-t GARCH input estimation.
   - Checks finite, symmetric covariance output and forecast-variance diagonals.
5. `test_cla_api`
   - Compiles and runs all original-name compatibility entry points.

The fixed three-asset references were independently reconstructed from the
published R recursion in NumPy. The expected turning points are:

```text
lambda          sigma          mu
0.416339869281  0.182756668825 0.102000000000
0.026683214985  0.080651550863 0.041291363591
0.023821645681  0.080437169956 0.040607580968
0.000000000000  0.080248818705 0.039337080713
```

The GNU Fortran strict configuration uses warnings as errors, bounds checking,
undefined-variable initialization checks where available, and floating-point
traps. The optimized configuration uses `-O3` and reruns every target.
