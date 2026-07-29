# API

All public names are available from module `evir`.

## Distribution functions

| Original R name | Fortran procedure | Purpose |
|---|---|---|
| `dgev` | `dgev` | GEV density, including the Gumbel limit |
| `pgev` | `pgev` | GEV CDF |
| `qgev` | `qgev` | GEV quantile |
| `rgev` | `rgev` | GEV random generation |
| `dgpd` | `dgpd` | GPD density, including the exponential limit |
| `pgpd` | `pgpd` | GPD CDF |
| `qgpd` | `qgpd` | GPD quantile |
| `rgpd` | `rgpd` | GPD random generation |

Scalar distribution functions are elemental and therefore also accept arrays.
Random generation uses an explicit `evir_rng` state.

## Model fitting

| Original R name | Fortran procedure | Result type |
|---|---|---|
| `gev` | `gev` | `gev_fit_result` |
| `gumbel` | `gumbel` | `gev_fit_result` with `gumbel=.true.` |
| `gpd` | `gpd` | `gpd_fit_result` |
| `pot` | `pot` | `pot_fit_result` |
| `gpdbiv` | `gpdbiv` | `gpdbiv_fit_result` |

`gev` and `gumbel` accept an optional integer `block_size`. `gpd`, `pot`, and
each bivariate margin accept exactly one of a numeric threshold or an integer
number of upper extremes. `gpd` supports `method='ml'` or `method='pwm'` and
`information='observed'` or `information='expected'`.

## Return levels and tail risk

| Original R name | Fortran procedure |
|---|---|
| `rlevel.gev` | `rlevel_gev`, `rlevel_gev_profile` |
| `gpd.q` | `gpd_q`, `gpd_q_wald`, `gpd_q_profile` |
| `gpd.sfall` | `gpd_sfall`, `gpd_sfall_profile` |
| `riskmeasures` | `riskmeasures` |

The profile procedures return `profile_result`, including the profile grid and
log-likelihood values.

## Exploratory and diagnostic algorithms

| Original R name | Fortran procedure | Result type |
|---|---|---|
| `emplot` | `emplot` | `xy_result` |
| `meplot` | `meplot` | `xy_result` |
| `qplot` | `qplot` | `xy_result` |
| `records` | `records` | `records_result` |
| `hill` | `hill` | `band_result` |
| `exindex` | `exindex` | `matrix_result`, columns N, K, threshold, theta2, theta |
| `shape` | `shape` | `band_result` |
| `quant` | `quant` | `band_result` |
| `tailplot` | `tailplot` | `tail_curve_result` |
| `decluster` | `decluster` | `decluster_result` |
| `findthresh` | `findthresh` | scalar threshold |

## Bivariate interpretation

- `interpret_gpdbiv` returns marginal exceedance probabilities, joint
  exceedance probability, independence product, and both conditional
  exceedance probabilities.
- `bivariate_cdf` evaluates the fitted logistic bivariate CDF.
- `bivariate_survivor` evaluates the joint survivor probability.
- `logistic_exponent` evaluates the exponent measure.

## Status codes

- `evir_ok`
- `evir_invalid_input`
- `evir_no_exceedances`
- `evir_optimization_failed`
- `evir_singular_hessian`
- `evir_domain_error`
