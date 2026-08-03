# Validation

The permanent test suite exercises:

1. PCA orthogonality, spanned-factor construction, and synthetic VAR recovery.
2. Foreign/star factors, GVAR assembly, and transition-weight normalization.
3. Closed-form scalar affine recursions, factor rotations, Gaussian likelihood,
   and measurement-covariance mapping.
4. Forecasts, orthogonal/generalized responses, FEVD/GFEVD normalization,
   expected short rates, term premia, and forward rates.
5. BFGS/Nelder-Mead optima, numerical derivatives, stationarity projection,
   and PSD/block covariance transforms.
6. Reproducible IID and block resampling, percentile bounds, VAR bootstrap,
   stochastic-approximation bias correction, and stability shrinkage.
7. JLL factor reconstruction, restriction masks, and covariance symmetry.

The final archive is built in two modes:

```text
checked:   -O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
optimized: -O3 -std=f2018 -Wall -Wextra -Werror
```

Both modes link against BLAS and LAPACK and run all tests, examples, and the
demonstration executable.
