# API

All procedures are provided by `use mixedindtests`. Real calculations use
`real(dp)`, where `dp = kind(1.0d0)`.

## Data preparation and pairwise dependence

- `preparedata(x)` returns sorted unique values, empirical CDF, and masses.
- `stat_dep(x, y)` returns multilinear Kendall tau, Spearman rho, and scale.
- `stat_dep_ser(x, lag)` applies the same statistics to a circular lag.
- `EstDep(x)` computes pairwise matrices and chi-square combinations.
- `EstDepSerial(x, lag)` computes lag 1 through `lag` statistics.

## Cramer-von Mises kernels

- `Sn_A(x, trunc_level)` handles an observation matrix.
- `Sn_Aserial(x, p, trunc_level)` handles `p` consecutive scalar observations.
- `Sn_AserialVec(x, p, trunc_level)` handles consecutive random vectors.
- `Sn_serial(x, p)` returns only the global serial statistic.
- `bootstrap(multiplier, sn_multiplier, xi)` evaluates one centered multiplier
  realization.

The `sn_result` type contains `stats`, `cardinality`, binary `subsets`, the
three-dimensional multiplier array, the global `sn`, and its multiplier
matrix.

## Bootstrap tests

- `TestIndCopula(x, trunc_level, b, seed)` tests mutual independence.
- `TestIndSerCopula(x, p, trunc_level, b, seed)` tests scalar randomness.
- `TestIndSerCopulaMulti(x, p, trunc_level, b, seed)` tests vector randomness.

Optional arguments use upstream defaults (`trunc_level=2`, `b=1000`). Results
are percentages and are returned in `copula_test_result`.

## Moebius covariance tests

- `EstDepMoebius(x, trunc_level)`
- `EstDepSerialMoebius(y, p, trunc_level)`

The result contains normalized Spearman, van der Waerden, and Savage
coefficients, full and pair-only chi-square combinations, p-values,
cardinalities, and subset indicators.

## Selection and simulation

- `select_p(x, p0, d, q, lambda)` selects a serial order.
- `SimAR1Poisson(param, n, seed)` simulates the conditional-Poisson AR(1).
- `SimCopulaSeries(family, n, tau, param, seed)` supports `ind`, `tent`,
  `gaussian`, `t`, `clayton`, `fgm`, `frank`, `gumbel`, `joe`, and `plackett`.
  `param` is the Student degrees of freedom or the FGM Markov order.
- `Finv(u, k)` evaluates the seven upstream simulation margins.

## Status values

Result types contain `status`:

- `mixedind_success`
- `mixedind_invalid_argument`
- `mixedind_allocation_error`
- `mixedind_numerical_error`
