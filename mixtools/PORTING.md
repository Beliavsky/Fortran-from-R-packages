# Porting notes

## Interface mapping

R vectors, matrices, lists, formulas, and S3 objects are replaced by explicit
Fortran arrays and derived types. Fortran is case-insensitive, and R dots are
mapped to underscores.

IEEE NaNs are not used as a general missing-data interface in this package.
Callers should remove or impute nonfinite observations before fitting.

## Closely translated algorithms

The following preserve the principal upstream likelihood and EM/ECM update
structure:

- univariate normal mixtures
- multivariate normal mixtures
- gamma mixtures
- multinomial mixtures
- repeated-measures normal mixtures
- linear, logistic, and Poisson regression mixtures
- exponential and Weibull right-censored mixtures
- normal-mixture linear mean/inverse-variance constraints
- weighted kernel densities, weighted quantiles, and simulation helpers

## Deliberate adaptations

- `npMSL` uses the common deterministic product-kernel mixture engine rather
  than reproducing the full R bandwidth-optimization loop.
- `spRMM_SEM` uses a flexible censored Weibull basis as a deterministic
  semiparametric reliability approximation.
- `weibullRMM_SEM` uses posterior expectations rather than stochastic label
  draws. It fits the same censored Weibull mixture likelihood.
- `regmixEM.mixed` is represented by ridge-regularized group intercepts rather
  than the complete random-effect covariance machinery of the R routine.
- `segregmixEM` uses safeguarded hinge-design breakpoint exploration instead
  of calling the R `segmented` package.
- `flaremixEM` is represented by a barrier-stabilized Gaussian regression
  mixture. The specialized upstream FLARE error model is not reproduced.
- `hmeEM` retains the upstream two-component limitation and uses logistic
  gating with weighted expert regressions.
- `regmixMH` is a self-contained Gibbs sampler with weak default priors rather
  than an exact reproduction of every upstream hyperparameter option.
- `boot.comp` and `boot.se` are implemented for normal mixtures. The original
  R dispatch across many model families is not emulated.
- equality tests use likelihood-ratio statistics with chi-square calibration;
  the R bootstrap path remains available separately for normal mixtures.

## External R dependencies

- `kernlab`: plotting/dimension-reduction use is omitted.
- `MASS`: required linear algebra and random generation are implemented in
  Fortran.
- `segmented`: replaced by the typed hinge-design search noted above.
- `survival`: censored likelihoods are implemented directly.
- `plotly` and `scales`: display-only and omitted.

## Numerical safeguards

All covariance solves use Cholesky-based routines with diagonal stabilization.
Scales are bounded away from zero. Posterior probabilities are normalized in
log space. Invalid arguments and singular systems return status codes rather
than terminating the process.
