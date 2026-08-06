# Validation

Six deterministic test programs cover:

1. Kernel values, convolution, boundary compatibility, eigenanalysis, and linear solves.
2. Local PCA dimensions, nonnegative weights, factor-selection IC values, and reconstruction monotonicity.
3. Adaptive POET rho selection, diagonal preservation, symmetry, and positive definiteness.
4. Analytical portfolio identities, ARMA expected-return forecasts, and complete prediction workflow.
5. Tapered bootstrap covariance and the Su-Wang statistic with seeded bootstrap replication.
6. Expanding-window weight budgets, holding periods, return lengths, and annualized metrics.

The demonstration estimates portfolio weights and compares rolling TV-MVP and equal-weight Sharpe ratios on deterministic synthetic returns.

Validated compiler configuration:

- GNU Fortran 14.2.0
- Checked: `-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999`
- Optimized: same language/warning policy with `-O3`
