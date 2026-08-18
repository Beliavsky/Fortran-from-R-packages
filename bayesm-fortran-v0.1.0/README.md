# bayesm-fortran

Modern free-format Fortran translation of the computational core of the R
package `bayesm` 3.1-7.

The library provides Bayesian regression/SUR, binary/ordinal/multinomial and
multivariate probit, MNL and nonhomothetic-logit utilities, negative-binomial
samplers, finite and Dirichlet-process normal mixtures, hierarchical linear
and discrete-choice samplers, instrumental-variable Gibbs sampling,
scale-usage latent-variable sampling, and aggregate random-coefficient logit
(BLP) routines.

## Build

```text
fpm build
fpm test
fpm run --example bayesm_example
```

All compiled source is free-format `.f90`.  The library is self-contained and
has no R, Rcpp, Armadillo, BLAS/LAPACK, C, or C++ runtime dependency.

## API naming

R names are mapped to Fortran snake_case names, for example:

- `createX` -> `create_x`
- `lndMvn` -> `lnd_mvn`
- `rmnpGibbs` -> `rmnp_gibbs`
- `rhierLinearModel` -> `rhier_linear_model`
- `rhierMnlRwMixture` -> `rhier_mnl_rw_mixture`
- `rscaleUsage` -> `rscale_usage`
- `rbayesBLP` -> `rbayes_blp`

The top-level `bayesm` module re-exports the public computational API and its
result/data derived types.

See `TRANSLATION_NOTES.md` for the complete coverage map and differences from
R list/formula/S3 interfaces.
