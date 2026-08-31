# API coverage

Upstream: R package `pan` version 2.0, published 2026-06-30.

## Exported R functions

| Upstream API | Fortran API | Coverage |
| --- | --- | --- |
| `pan()` | `pan_mcmc` | Implemented numerical model and Gibbs updates |
| `pan.bd()` | `pan_bd_mcmc` | Implemented block-diagonal random-effect model |
| `ecme()` | `ecme_fit` | Implemented ML target, GLS special case, EB effects |

## `pan()` parity

Implemented:

* multivariate response;
* arbitrary fixed-effect design selected from `pred`;
* arbitrary random-effect design selected from `pred`;
* full `q*r` by `q*r` covariance for response-major vectorized random effects;
* `a`, `Binv`, `c`, and `Dinv` inverse-Wishart prior interpretation;
* initialization when no restart state is supplied;
* restart from beta/sigma/psi/completed-y state;
* partially missing rows;
* completely missing rows;
* chain storage for beta, sigma, and psi;
* final imputed response and final restart state.

Numerical restructuring relative to legacy `pan.f`:

* the 66 fixed-form workspace routines are collapsed into typed kernels;
* dense SPD operations use direct lower-Cholesky code rather than the legacy
  upper-triangular workspace convention;
* the package-local RNG is deterministic Park-Miller plus double-precision
  Box-Muller/Marsaglia-Tsang generation, rather than reproducing the legacy
  single-precision gamma implementation bit-for-bit.

These changes preserve the statistical distributions and conditioning
structure but do not promise bitwise identity to the historical Fortran RNG
stream.

## `pan.bd()` parity

Implemented:

* response-specific `c(j)` degrees of freedom;
* response-specific `Dinv(:,:,j)` scales;
* independent inverse-Wishart random-effect covariance blocks;
* the same missing-data, fixed-effect, residual-covariance, chain, and restart
  behavior as `pan_mcmc`.

## `ecme()` parity

Implemented:

* univariate response;
* fixed and random effects selected from `pred`;
* identity residual occasion covariance by default;
* extraction of cluster `V_i` from user-supplied `vmax` and `occ`;
* exact GLS when no random effects are requested;
* Gaussian mixed-model ML fitting when random effects are present;
* log-likelihood trace;
* convergence flag and iteration count;
* covariance of fixed effects;
* empirical Bayes random-effect means and conditional covariances.

Difference:

The upstream fixed-form routine `ecme3` uses Schafer's ECME sequence. The
modern implementation uses an EM/ECME-target iteration for the same Gaussian
marginal maximum-likelihood objective. It therefore targets the same ML
estimator but does not reproduce every intermediate upstream iterate.

## Intentionally omitted R-specific code

* `.Fortran` argument packing and workspace arrays;
* S3/list construction and R storage-mode coercion;
* `cat()` progress output;
* example plotting and `acf`;
* R datasets and `.rda` serialization;
* native-routine registration skeleton;
* R help/vignette build machinery.

These omissions do not remove a numerical model exported by the upstream
package.
