# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018
- Linux x86-64
- No external numerical library required

## Debug configuration

```text
-std=f2018 -ffree-line-length-none -Wall -Wextra
-Wimplicit-interface -Werror -fbacktrace -O0 -g
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

## Release configuration

```text
-std=f2018 -ffree-line-length-none -Wall -Wextra
-Wimplicit-interface -Werror -fbacktrace -O2
```

## Executed suites

1. Return, portfolio, and drawdown identities
2. Risk and performance-ratio calculations
3. CAPM, timing regression, co-moment, and portfolio-risk calculations
4. Rolling and cleaning calculations
5. Structured/plug-in M2/M3/M4, EWMA higher moments, MCA, NCE, GPD, Monte Carlo, lognormal/kernel risk, bootstrap errors, dynamic CAPM, capture, outperformance, and summary calculations
6. Exact finite-sample M2/M3/M4 shrinkage corrections, target covariances, unbiased coskewness MSE, and constrained multi-target workflows

## Exact finite-sample correction tests

- M2, biased M3, and M4 sample-loss and target-covariance terms are checked against independently coded empirical influence-function calculations.
- The unbiased M3 k-statistic path is checked against a fixed sixth-order reference vector.
- Every M2, M3, and M4 target family is exercised, including two observed factors.
- Deterministic b-vector references cover independent, equal-marginal, factor, constant-correlation, Simaan, and central-symmetry targets.
- Single-target weights are checked against the closed-form clipped b/A solution.
- Multi-target weights are checked for nonnegativity, sum at most one, symmetric A matrices, convergence, and finite estimates.
- Bias-corrected independent structured coskewness is checked against its exact n^2/((n-1)(n-2)) multiplier.

## Advanced test details

- Factor-generated three-asset data exercise observed and latent moment structures.
- Single- and multi-target shrinkage weights are checked for nonnegativity and simplex feasibility.
- MCA directions are checked for unit norm and nonzero higher-moment reconstruction.
- NCE covariance symmetry, positive residual variances, and skewness/kurtosis feasibility inequalities are checked.
- A deterministic GPD quantile grid checks shape and scale recovery and VaR/ES ordering.
- Monte Carlo VaR/ES uses reproducible seeds; portfolio ES contributions are checked for additivity.
- Moving-block bootstrap standard errors are checked for finite nonnegative output.
- Rolling, expanding, and conditional CAPM are checked against exact synthetic coefficients.
- Kernel VaR/ES contributions, lognormal risk, capture ratios, and trailing outperformance probabilities are exercised.

## Result

Both debug and optimized builds passed every suite and application. The GPL-2.0-or-later source-header audit also passed.
