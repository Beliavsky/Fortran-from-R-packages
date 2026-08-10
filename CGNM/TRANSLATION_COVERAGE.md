# Translation coverage

## Directly translated numerical behavior

The following comes from `R/Cluster_Gauss_Newton_method.R` and is represented
in the Fortran implementation:

- cluster initialization and invalid-model resampling;
- per-member target/weight matrices;
- algorithm version 3 local clustered linear approximation;
- algorithm version 1 global-distance approximation;
- reciprocal standardized-distance weights raised to `gamma`;
- regularized CGNR step solve;
- accept-if-SSR-decreases rule;
- lambda `/10` on acceptance and `*10` on rejection;
- `10^10 * initial_lambda` inactive-member cutoff;
- initial-range step restriction;
- parameter bound transformations;
- multi-objective pseudo-observations;
- bootstrap types 1-3;
- EBE observation weighting;
- core accepted/best-index postprocessing.

## Standalone replacements for external R functionality

- `stats::kmeans` -> Lloyd k-means.
- `MASS::ginv` -> symmetric Jacobi eigendecomposition pseudoinverse.
- R RNG -> intrinsic Fortran RNG. Seeds are reproducible within the Fortran
  implementation but do not reproduce R's random stream.
- `stats::qt` in Grubbs screening -> standalone inverse-normal plus finite-df
  expansion.

These substitutions can change exact trajectories while preserving the
translated CGNM equations.

## Intentional correction

The current R wrapper computes the inverse lower-bound transformation as
`log(theta-lower)` and reports the forward transformation as
`exp(x)+lower`, but one internal model wrapper contains `exp(x)-lower`.
Those three expressions cannot all be mutually inverse. The Fortran version
uses the mathematically consistent pair

```text
z = log(theta - lower)
theta = lower + exp(z)
```

and documents the difference instead of propagating a bound-violating sign
error.

## Not translated

The following are R presentation/orchestration layers rather than the core
numerical CGNM iteration and are intentionally omitted:

- Shiny application code;
- ggplot2 plotting functions;
- R formula/string `eval(parse())` reparameterization;
- `.RDATA` logging and profile-likelihood routines that reconstruct datasets
  by reading those R log files;
- S3/data-frame formatting and print helpers;
- k-means visualization/diagnostic objects.

Postprocessing functions whose sole purpose is to build plotting data are not
ported. Their underlying final-cluster, residual, target and history arrays are
available in `cgnm_result` for native Fortran analysis.
