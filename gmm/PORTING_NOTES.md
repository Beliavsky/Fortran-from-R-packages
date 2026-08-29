# Porting notes

## Translation policy

Only computational/statistical code was translated. R formula parsing, model frames, S3 object presentation, plotting, printing, and similar UI glue were omitted. Numerical routines are exposed through explicit arrays, callbacks, and derived result types.

The package's native `lambda_met.f` routines and R-level estimator mathematics were used as the algorithmic source. The upstream computational sources used for the translation are retained under `upstream/`.

## Supplied R helper module

The project reuses the supplied MIT-licensed `r_mod.f90` for applicable R-compatible probability, RNG, optimization, and time-series helpers. `upstream/r_mod-original.f90` is the exact supplied file. `src/r_mod.F90` differs only by free-form continuation/line wrapping needed by the build. No statistical helper from `r_mod` was reimplemented merely for convenience.

Project-specific helpers in `gmm_linalg` and related modules cover genuinely missing GMM/GEL matrix algebra and estimating-equation operations.

## GMM

`gmm_fit` accepts pure callback procedures returning the observation-by-moment matrix and the Jacobian of mean moments. Two-step, iterative, and CUE estimates are implemented. Weight matrices support MDS, iid, HAC, and identity paths. `gmm_fit_fixed_weight` covers the upstream fixed-W computational case.

Because the supplied `r_mod::optim_bfgs` objective interface is a pure one-argument function, generic GMM/GEL/ATE-GEL optimizers use module-local state to bridge user callbacks into the optimizer. Consequently those outer optimization calls are not reentrant or thread-safe while an optimization is active. Independent fits executed serially are unaffected.

## GEL

EL, ET, CUE, HD, ETEL, ETHD, and RCUE are implemented. The original Wu/Newton-style feasibility handling and CUE/active-set ideas are represented directly in modern Fortran. The outer GEL estimator computes implied probabilities, parameter covariance, lambda covariance, and J/LM/LR-style specification statistics.

The Fortran implementation exposes the numerical smoothing/covariance quantities directly. It does not reproduce R `tskernel` or `sandwich` object semantics.

## HAC and bandwidths

The covariance module contains the kernel and HAC subset needed by `gmm`; it is not a translation of the entire R `sandwich` package. Andrews bandwidth selection supports AR(1) and ARMA(1,1) approximations using the existing `r_mod` time-series fits.

The upstream `bwWilhelm` ARMA branch appears to compute `arma.coef` and subsequently refer to an `ar.coef` name. The Fortran port uses the fitted ARMA coefficients consistently; this is treated as an obvious variable-name defect rather than reproducing an undefined-name path.

## Linear and system estimators

Linear GMM/IV and 2SLS are direct matrix estimators. System GMM includes MDS/HAC/conditionally homoskedastic weighting and SUR, 3SLS, FIVE-style, and random-effects wrappers.

The numerical system API uses rectangular arrays `x(n,k,neq)` and `z(n,q,neq)`, so all equations in one call currently share the same numbers of regressors and instruments. Upstream R list/formula interfaces can represent equation-specific dimensions and elaborate cross-equation coefficient restrictions; that constructor/object layer is not recreated here. Common versus equation-specific coefficient vectors are supported by the Fortran numerical core.

## ATE GEL

`ate_moments` and `ate_gradient` reproduce the estimating-equation kernels with a packed numeric data matrix. `ategel_fit` connects them to the generic GEL estimator. Linear/logit/probit outcome links and balance/sample-balance/ATT moment variants are provided. R formula/family objects are not required.

## Tests and validation

The retained tests cover:
- independent-reference linear two-step, iterated, CUE GMM and 2SLS coefficients;
- a generic nonlinear callback GMM problem;
- EL, ET, HD, CUE, ETEL, and ETHD estimates compared with independent NumPy/SciPy nested optimizations;
- conditionally homoskedastic system GMM against an independent matrix calculation;
- stable characteristic-function reference values;
- GEL lambda/probability identities and HAC/kernel calculations;
- ATE moment/Jacobian formulas;
- Kleibergen K/J/S block formulas.

FPM was not available in the validation environment, so the same FPM source graph was compiled directly with GNU Fortran, BLAS, and LAPACK using Fortran 2018, runtime checks, and implicit-interface errors.
