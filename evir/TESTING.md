# Testing

The release is tested in two configurations with GNU Fortran 14.2:

## Strict build

```text
-std=f2018
-Wall -Wextra -Werror
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-fbacktrace
-O0 -g
```

## Optimized build

```text
-std=f2018
-Wall -Wextra -Werror
-O3
```

The four test programs cover:

1. GEV/GPD density, CDF, quantile inversion, zero-shape limits, and RNG support.
2. Threshold selection, declustering, empirical tails, mean excess, QQ data,
   records, Hill estimates, and extremal-index tables.
3. GEV/Gumbel/GPD/PWM/point-process fitting, expected covariance, return levels,
   profile intervals, Wald intervals, VaR, and expected shortfall.
4. Local bivariate fitting, logistic CDF and survivor probabilities, and
   conditional exceedance interpretation.

Independent Python/SciPy references include:

- deterministic GEV MLE `(xi, sigma, mu) =`
  `(0.1497925024, 1.99375116, 5.00003276)`;
- deterministic GPD MLE `(xi, beta) =`
  `(0.23939195, 1.20947820)`;
- fixed GEV, GPD, and point-process negative log likelihoods;
- fixed bivariate logistic CDF and survivor probabilities.
