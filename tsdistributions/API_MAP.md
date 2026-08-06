# API map

| R API | Fortran API | Status |
|---|---|---|
| `ddist`, `pdist`, `qdist`, `rdist` | same names | Implemented |
| `d/p/q/rnorm` through wrapper | `ddist/pdist/qdist/rdist('norm',...)` | Implemented |
| `d/p/q/rstd` | same names | Implemented |
| `d/p/q/rsnorm` | same names | Implemented |
| `d/p/q/rsstd` | same names | Implemented |
| `d/p/q/rged` | same names | Implemented |
| `d/p/q/rsged` | same names | Implemented |
| `d/p/q/rnig` | same names | Implemented |
| `d/p/q/rgh` | same names | Implemented |
| `d/p/q/rjsu` | same names | Implemented |
| `d/p/q/rghst` | same names | Implemented |
| `d/p/q/rghyp` | same names plus `*_raw` | Implemented |
| `nigtransform`, `ghyptransform` | same names | Implemented |
| `dskewness`, `dkurtosis` | same names | Implemented |
| `authorized_domain` | same name | Numerical grid equivalent |
| `distribution_modelspec` | same name | Typed specification |
| `distribution_bounds` | same name | Implemented |
| `estimate.tsdistribution.spec` | `estimate_distribution` | Numerical MLE equivalent |
| `coef`, `logLik`, `AIC`, `BIC` | fields of `distribution_fit` | Implemented |
| `bread`, `estfun`, `vcov` | score/Hessian fields and `information_covariance` | Implemented |
| `tsmoments` | `distribution_moments` | Implemented |
| `tsprofile` | same name | Implemented |
| `spd_modelspec` | same name | Typed specification |
| `estimate.tsdistribution.spdspec` | `estimate_spd` | PWM/KDE equivalent |
| `dspd`, `pspd`, `qspd`, `rspd` | same names | Implemented |
| plot/print/summary methods | Documentation/client code | Omitted |
