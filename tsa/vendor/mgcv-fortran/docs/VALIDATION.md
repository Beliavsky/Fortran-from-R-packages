# Validation

The release test program covers:

- exact inverse relationships for `notExp`/`notLog` and
  `notExp2`/`notLog2`;
- SPD inversion, symmetric matrix roots, and tridiagonal Cholesky;
- cubic, P-spline, cyclic, prediction, and tensor-product basis dimensions;
- basis centering and reconstruction on training points;
- Gaussian P-spline GAM fitting with automatically selected GCV smoothing;
- a fixed-lambda ridge result independently computed as
  `beta = solve(X'X + lambda S, X'y) = (2, 1)`;
- Poisson and binomial penalized IRLS;
- monotone penalized constrained least squares;
- multivariate-normal simulation moments;
- the 0.95 quantile of a one-degree-of-freedom chi-square distribution via
  `weighted_chisq_cdf`;
- nonnegative compound Poisson-gamma Tweedie draws; and
- deterministic `gamSim` dimensions and noise scale.

The synthetic Gaussian GAM test uses

```text
f(x) = sin(2*pi*x) + 0.4*cos(6*pi*x)
```

with deterministic high-frequency contamination. The required fitted RMSE is
below 0.18 and the selected effective degrees of freedom must lie strictly
between the parametric minimum and full basis dimension.

Build configurations used for the release:

```text
Debug:     -std=f2018 -Wall -Wextra -Wimplicit-interface
           -fcheck=all -fbacktrace -O0 -g
Optimized: -std=f2018 -O3 -march=native
```

The final archive is extracted into a fresh directory and both the test suite
and example are rebuilt from that extracted copy.
