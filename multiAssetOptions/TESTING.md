# Testing

The test suite contains six independent programs:

1. `test_grid_payoff` checks uniform/nonuniform spacing, mesh shifting, and all
   payoff families.
2. `test_operator` compares a one-dimensional operator against hand-derived
   stencil coefficients and checks the discounted-constant identity.
3. `test_european` compares call and put values with Black-Scholes formulas and
   checks put-call parity.
4. `test_american` checks the early-exercise premium, an accepted benchmark
   range, the intrinsic-value floor, and actual penalty iteration.
5. `test_multiasset` compares an independent two-asset digital option against
   the product of analytic marginal probabilities.
6. `test_adaptive` checks adaptive-step completion, ordering, and pricing
   accuracy.

Run with FPM:

```sh
fpm test
```

Or use the scripts in `scripts/` for checked and optimized GNU Fortran builds.
