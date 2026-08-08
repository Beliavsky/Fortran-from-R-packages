# Reference generation

The validation suite does not rely on the translated implementation to generate its own expected answers.

## Analytical references

- Black-Scholes-Merton values and Greeks use closed-form expressions.
- Bond prices, duration, convexity, and yield repricing use direct discounted-cash-flow identities.
- Portfolio tests check budget, box, group, target-return, and risk-contribution identities directly.
- `random_returns(exact=.true.)` is checked against requested finite-sample means, standard deviations, and correlation matrices.
- Marginal risk contributions are compared with central finite differences.
- CPPI is checked through its cash/bond/exposure accounting identities.

## Original package references

The retained original tree contains:

- `original/inst/tinytest/`
- `original/inst/unitTests/`
- package examples and manual pages

Representative fixed values and invariants from these materials were transcribed into the Fortran tests where they test a computational routine rather than R presentation behavior.

## Independent numerical references

- Gauss-Legendre and Gauss-Chebyshev nodes and weights are checked against defining moment identities.
- Heston and barrier prices are checked against published/reference values used by the original package tests.
- Optimization methods are run on benchmark functions with known or easily independently verified minima.
- Matrix repair is checked through symmetry and nonnegative eigenvalues.
- PBO is checked through explicit combinatorial splits and rank/logit construction.
