# bivpois-fortran

Modern Fortran translation of the computational code in the R package **bivpois** 1.2 by Michail Tsagris.

The project is organized as an FPM library and requires no external numerical libraries. Plotting and R-specific presentation code are not translated.

## Implemented functionality

- Bivariate Poisson log-PMF and PMF (`dbp`, `dbp_scalar`).
- Random generation by the three-independent-Poisson construction (`rbp`).
- Fast profile-likelihood MLE corresponding to R `bp.mle2` (`bp_mle2`).
- Full MLE/inference corresponding to R `bp.mle` (`bp_mle`), including:
  - three lambda estimates,
  - implied correlation,
  - independence and fitted log-likelihoods,
  - likelihood-ratio p-value,
  - observed-information and asymptotic Wald calculations,
  - 95% confidence intervals.
- Numerical profile likelihood and profile-grid confidence interval (`lambda3_profile`).
- Loukas-Kemp/Rayner-Thas-Best dispersion statistic and parametric-bootstrap GOF test (`bp_gof`, `bp_gof2`).
- Contingency-table construction (`make_bp_table`).
- Probability-grid computation replacing the computational part of `bp.contour` (`bp_probability_grid`).
- Reproducible RNG seeding (`seed_rng`).

## Numerical changes

The R density code evaluates the Karlis-Ntzoufras formula through binomial coefficients and powers of `lambda3/(lambda1*lambda2)`. The Fortran port evaluates the equivalent latent-Poisson convolution in log space using log-sum-exp. This is more robust for large counts and remains defined when one or more component rates are zero.

The MLE remains the same one-dimensional profile problem: the marginal means determine `lambda1 + lambda3` and `lambda2 + lambda3`, leaving only `lambda3` to optimize. The R package's `min(mean1,mean2)-0.05` upper-margin convention is preserved when feasible. Sparse cases for which that interval would be empty fall back to the mathematically valid nonnegative interval.

R's `bp.gof2` exists to remove R interpreter overhead by vectorizing Monte Carlo replicates. In compiled Fortran, the loop itself is already efficient, so `bp_gof2` is an API-compatible alias of the same memory-efficient implementation as `bp_gof`.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example demo_bivpois
```

The validation environment did not contain FPM, so the exact source/test units were also compiled directly with GNU Fortran 14.2 using Fortran 2018 mode, warnings, implicit-interface errors, and runtime checking. See `VALIDATION.md`.

## Modules

- `bivpois`: umbrella public API.
- `bivpois_distribution`: PMF, simulation, probability grids, contingency tables.
- `bivpois_fit`: likelihood, MLE, inference, profile likelihood.
- `bivpois_gof`: dispersion statistic and Monte Carlo GOF.
- `bivpois_math`: RNG and small statistical utilities.
- `bivpois_kinds`: floating-point kind definition.

## Licensing and provenance

The upstream R package declares `GPL (>= 2)`. This translation is therefore distributed under **GPL-2.0-or-later**. The original uploaded package archive and source/metadata snapshots are retained under `provenance/`.
