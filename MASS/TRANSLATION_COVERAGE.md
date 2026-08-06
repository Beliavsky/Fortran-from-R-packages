# Translation coverage

## Implemented

The computational kernels behind `area`, `bandwidth.nrd`, `bcv`, `boxcox`,
`con2tr`, `contr.sdif`, `corresp`, `cov.trob`, `dose.p`, `fbeta`, `fitdistr`,
`fractions`, `gamma.shape`, `ginv`, `glm.nb`, `huber`, `hubers`, `isoMDS`,
`kde2d`, `lda`, `lm.gls`, `lm.ridge`, `loglm`, `logtrans`, `mca`, `mvrnorm`,
`nclass.freq`, negative-binomial family calculations, `negexp.SSival`, `Null`,
`polr`, psi functions, `qda`, `rational`, `rlm`, `rnegbin`, `sammon`, `Shepard`,
residual diagnostics, theta estimators, `ucv`, `width.SJ`, robust covariance,
and LQS/LMS/LTS regression are represented.

`addterm`, `dropterm`, and `stepAIC` have explicit design-matrix equivalents.
Histogram bin-count and plotting wrappers are not duplicated; the bin-count
rules needed by other routines are available numerically.

## Omitted by design

- Plotting: `eqscplot`, `frequency.polygon`, `ldahist`, `parcoord`, `truehist`,
  and all plot/pairs/biplot methods.
- Formula, S3, printing, summary, profile, update, and model-frame machinery.
- Dataset objects and examples that only demonstrate R interfaces.
- `glmmPQL`, which requires `nlme` mixed-effects model classes and methods.
- `rms.curv`, whose input is an R nonlinear-model object carrying gradient and
  Hessian attributes.
- Formula renumbering/list utilities (`denumerate`, `renumerate`, `enlist`).

These omissions do not remove kernels used by the implemented array APIs.
