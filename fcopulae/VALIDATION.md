# Validation

## Environment

- GNU Fortran: 14.2.0
- Language mode: Fortran 2018
- Libraries: LAPACK and BLAS

Debug flags:

```text
-O0 -g -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

Release flags:

```text
-O2 -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

## Test suites

### `test_archimedean`

- Default parameter validity for all 22 families
- Generator/inverse identity for all 22
- Frechet CDF bounds and finite nonnegative densities
- Random generation in the unit square for all 22
- Independence and Clayton closed-form identities
- Clayton and Gumbel Kendall-tau references
- Clayton Spearman-rho and tail-dependence checks
- Executed bounded fitting path for every family

### `test_elliptical`

- Bivariate Normal quadrant-probability identity
- Kendall and Normal Spearman formulas
- Marginal quantile/CDF inversion for all seven families
- Positive marginal densities and finite copula densities
- CDF bounds and random generation for all seven
- Approximate uniformity of simulated margins
- Normal and Student-t tail-dependence checks
- Executed bounded fitting path for every family

### `test_extreme_empirical`

- Pickands bounds for all five official families
- Finite CDF/density, simulation, tau, rho, and tail paths for all five
- Gumbel equivalence between extreme-value and Archimedean representations
- Gumbel analytical tau and upper-tail coefficient
- Executed bounded fitting path for every official family
- Supplemental Gumbel-II, independence, and comonotonic Pickands evaluations
- Empirical copula corners, point probability, and density mass
- Frechet and Marshall-Olkin references
- Debye `D_1(1)` reference value

## Applications exercised

- `demo_fcopulae`
- `dependence_example`
- `fit_csv` with Archimedean type 1
- `fit_csv` with the Gaussian elliptical copula
- `fit_csv` with the Gumbel extreme-value copula

## Interpretation limits

The tests establish executable numerical functionality and selected analytical identities. They do not claim iteration-for-iteration equivalence with R's optimizers or random-number generators. Nonstandard elliptical simulation uses rank transforms and is tested for copula-range and marginal-uniformity behavior, not exact R draw equivalence.
