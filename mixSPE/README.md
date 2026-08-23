# mixSPE-fortran

Modern Fortran/FPM translation of the numerical core of the R package `mixSPE`.

## Implemented API

- `dpe`, `log_dpe`: multivariate power-exponential density.
- `rpe`: multivariate power-exponential simulation.
- `cov_pe`: covariance corresponding to a power-exponential scale matrix.
- `dspe`, `log_dspe`: skew power-exponential density.
- `rspe`: Metropolis-Hastings sampler matching the upstream proposal design and
  using the supplied `mvtnorm-fortran` normal generator/density.
- `em_fit`: generalized EM fitting for a requested mixSPE model.
- `emgr_fit`: fit a grid of group counts and model names and select by BIC.
- `spe_model`: fitted-model type containing mixing proportions, locations,
  eigenvalues/eigenvectors, shape parameters, skew parameters, posterior
  memberships, MAP labels, log likelihood, BIC, and iteration count.
- `model_num_parameters`, `map_labels`.

Supported covariance model prefixes are `EII`, `VII`, `EEI`, `VVI`, `EEE`,
`EEV`, `VVE`, and `VVV`. The fourth character controls beta (`E`, `V`, or `D`),
and a fifth character enables skewness, following the upstream package.

## Build

```sh
fpm test
```

or compile with a Fortran 2018 compiler.

## Numerical translation notes

The density normalization and beta Newton update are direct translations of the
R formulas. The covariance structures are enforced through eigendecomposition
of weighted component scatter matrices. For the coupled eigen-structure models
(`EEV` and `VVE`) the Fortran implementation uses a deterministic generalized-EM
projection onto the requested constraint family instead of reproducing the R
package's Stiefel-manifold Armijo optimizer line-for-line. This keeps the same
model family and yields a valid generalized EM step, but fitted values need not
be bit-for-bit identical to R.

Initialization uses deterministic farthest-point/Lloyd k-means rather than R's
random `kmeans`/exponential-weight multi-start machinery. `emgr_fit` can be
called repeatedly with user-selected initial posterior weights via `em_fit` if
multi-start behavior is desired.
