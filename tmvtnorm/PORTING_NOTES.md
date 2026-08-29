# Porting notes

## Scope

The port targets the computational/statistical content of tmvtnorm 1.7. Plotting, examples embedded in R files, S3/S4 presentation, R `methods`/`Matrix` class dispatch, `stats4::mle` objects, and native registration/C wrappers are omitted.

The original computational R and Fortran files are retained in `upstream/` for provenance and algorithm comparison.

## mvtnorm dependency

The user-supplied `mvtnorm-fortran` port is included unchanged as `vendor/mvtnorm-fortran` and referenced by FPM as a path dependency. Its probability controls and numerical implementations are therefore the probability engine for this translation.

## Marginal densities

The upstream R implementation of `dtmvnorm.marginal2` derives an `(n-2)`-dimensional conditional rectangle using partial correlations. The Fortran implementation computes the same object directly with standard Gaussian conditioning:

`X_o | X_s=x_s ~ N(mu_o + Sigma_os Sigma_ss^-1 (x_s-mu_s), Sigma_oo-Sigma_os Sigma_ss^-1 Sigma_so)`.

The marginal density is then the untruncated one- or two-dimensional Gaussian density times this conditional rectangle probability divided by the total truncation probability. This is algebraically equivalent and avoids duplicated correlation-matrix manipulation.

## Gibbs sampling

For dense covariance input the upstream code explicitly forms each coordinate's conditional regression separately. The Fortran implementation inverts the covariance once and then uses

`E[X_i|X_-i] = mu_i - H_ii^-1 H_i,-i (X_-i-mu_-i)` and `Var[X_i|X_-i] = H_ii^-1`.

The same kernel is used for explicit precision input. General linear constraints compute each coordinate's feasible interval from every row of `D` before drawing from the corresponding univariate truncated conditional.

Sparse CSC/triplet samplers preserve the precision-form algorithm without constructing a dense matrix. CSC indices are exposed in native 1-based Fortran form.

## Student t sampler

`rtmvt_gibbs` follows the upstream Geweke construction: alternate a chi-square scale draw conditional on the current Gaussian latent state, truncated-normal coordinate updates, and the Student transformation `x = mu + z/w`.

## Moments

`mtmvnorm` ports the Manjunath-Wilhelm / Tallis / Leppard-Tallis first- and second-moment formulas. When not every coordinate is truncated it applies the Johnson-Kotz block reduction, as upstream does.

## Estimation

The R package builds dynamic named functions for `stats4::mle`. Fortran instead returns `tmvnorm_fit_t` and uses a small local Nelder-Mead optimizer with an impure callback, because each likelihood evaluation can invoke randomized/quasi-Monte-Carlo multivariate probability code.

The optional Cholesky parameterization uses log diagonal elements. This is intentionally more robust than allowing an unconstrained optimizer to propose singular Cholesky diagonals while preserving the same covariance model.

`gmm_tmvnorm` translates both upstream moment-condition systems. It performs identity-weight first-stage estimation followed by an IID moment-covariance weighted second stage. The external R `gmm` package's full option surface (HAC kernels, S3 objects, etc.) is dependency functionality rather than tmvtnorm code and is not recreated here.

## Numerical correction

For `dtmvt(..., log=TRUE)`, the mathematical truncated log-density is `log f_t(x) - log P(truncation)`. The upstream R source divides the log-density by the probability inside `ifelse`, which appears to be an implementation error. The Fortran port uses the mathematically correct log normalization.
