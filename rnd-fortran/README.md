# RND for modern Fortran

This project translates the computational routines of the CRAN package `RND`
1.2 into self-contained modern Fortran organized as an FPM package.

`RND` extracts option-implied risk-neutral densities using a single lognormal
model, a generalized-beta distribution, mixtures of lognormal distributions,
an Edgeworth expansion, Shimko's quadratic implied-volatility method, and an
ad hoc three-component lognormal mixture.

## Build and test

```text
fpm build
fpm test
fpm run rnd_demo
fpm run --example basic_pricing
fpm run --example implied_smile
```

The project has no external dependencies. It uses standard Fortran intrinsic
functions and includes its own regularized incomplete-beta calculation,
linear solvers, bounded-by-penalty objectives, Nelder-Mead optimizer, and
finite-difference Hessian.

When FPM is unavailable, GNU Fortran validation can be run with:

```text
./tools/validate.sh
```

## Main modules

- `rnd_densities`: generalized-beta, mixture-lognormal, Edgeworth, and Shimko densities.
- `rnd_pricing`: BSM, generalized-beta, mixture, Edgeworth, Shimko, and ad hoc option prices.
- `rnd_objectives`: translated calibration objective functions.
- `rnd_fitting`: model extraction, implied-volatility inversion, rates, point estimates, and `moe`.
- `rnd_optimize`: self-contained Nelder-Mead optimization and numerical Hessians.
- `rnd`: convenience module re-exporting the public API.

R function names containing dots use underscores in Fortran. For example,
`price.bsm.option` becomes `price_bsm_option`, and `extract.rates` becomes
`extract_rates`.

## Scope

All exported numerical R functions are represented. R graphics, PDF creation,
CSV side effects, formula objects, model summaries, and serialized `.rda` data
are not compiled. The original data files remain under `original/`.

The R `MOE` workflow is available as `moe` and `fit_all_densities`. It returns
the fitted model structures and finite-difference point density rather than
creating plots and CSV files. Callers can use the pricing and density routines
to generate any desired output grid.

## Licensing

The original package declares `GPL (>= 2)`. This translation therefore uses
`GPL-2.0-or-later`. See `LICENSE`, `NOTICE`, and the retained original source.
