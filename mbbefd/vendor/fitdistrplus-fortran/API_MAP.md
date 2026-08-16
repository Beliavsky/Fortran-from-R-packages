# API map

| R API | Fortran API | Status |
|---|---|---|
| `fitdist` | `fitdist`, `fitdist_auto` | Implemented for MLE/MME/QME/MGE/MSE |
| `fitdistcens` | `fitdistcens`, `mledist_censored` | Implemented for censored MLE |
| `mledist` | `mledist` | Implemented, weighted and bounded |
| `mmedist` | `mmedist` | Implemented for raw moments |
| `qmedist` | `qmedist` | Implemented with type-7 quantiles |
| `mgedist` | `mgedist` | All eight upstream distances |
| `msedist` | `msedist` | All five phi-divergence families |
| `descdist` | `descdist` | Numerical summary implemented; graph omitted |
| `gofstat` | `gofstat` | KS, CvM, AD, chi-square, AIC, BIC |
| `detectbound` | `detectbound` | Uses explicit model parameter domains |
| `prefit` | `prefit` | Bounded preliminary fit from feasible parameters |
| `bootdist` | `bootdist` | Parametric or nonparametric serial bootstrap |
| `bootdistcens` | `bootdistcens` | Censored-row bootstrap |
| `CIcdfplot` | `cdf_bootstrap_band` | Confidence-band computation only |
| `quantile.fitdist` | `dist%quantile(prob, fit%estimate)` | Direct callback equivalent |
| `logLik`, `AIC`, `BIC`, `coef`, `vcov` | fields of `fit_result` | Implemented |
| `Surv2fitdistcens` | `type(censored_sample)` | Replaced by typed left/right arrays |
| internal Turnbull/NPMLE utilities | `turnbull_npmle` | Deterministic endpoint-support EM equivalent |
| `plotdist`, `plotdistcens` | none | Plotting omitted |
| `cdfcomp`, `denscomp`, `qqcomp`, `ppcomp` | direct callback evaluation | Plotting wrappers omitted |
| censored comparison plots | none | Plotting omitted |
| `llplot`, `llcurve`, `llsurface` | objective callbacks in fitting engine | Graphics omitted |
| S3 print/summary/density methods | typed result fields | R-only infrastructure omitted |
