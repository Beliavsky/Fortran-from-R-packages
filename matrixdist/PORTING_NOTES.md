# Porting notes

## Translation strategy

The R package combines S4 classes with a large Rcpp/RcppArmadillo numerical layer.  This port does not emulate S4.  Model parameters are passed as ordinary Fortran vectors/matrices/cubes or small derived types.  Closely related native entry points (RK/PADE/uniformization variants and repeated distribution wrappers) are consolidated around shared, tested kernels rather than exposed as ~130 thin wrappers.

The matrix exponential is the upstream Pade-6 scaling/squaring algorithm.  Van Loan block exponentials are reused for PH EM occupancy integrals.  A uniformization implementation is also retained for parity/testing.

## Important documented corrections

1. **Discrete simulation absorption.**  Upstream `rdphasetype`/`rMDPHstar` constructs cumulative probabilities from only the transient `S` columns but loops until a declared absorbing state is reached.  With no absorbing column that state cannot be sampled.  The Fortran implementation explicitly appends the exit probability `1-S*1`, so the chain can absorb.

2. **Positive reward time scaling.**  The no-zero-reward branch of upstream `linear_combination` multiplies rows by rewards, while its zero-reward branch and standard PH reward time change imply division of rates by rewards.  `tvr_ph` and `linear_combination` consistently divide generator rows by positive rewards.

3. **Matrix-GEV CDF orientation.**  Upstream `mgevcdf` returns `1-alpha*exp(S*t)*1` although the GEV time transform `t=g^{-1}(x)` is decreasing.  In the one-state case this conflicts with both `mgevden` and `rmatrixgev`.  The Fortran CDF uses `alpha*exp(S*t)*1`, producing the standard GEV CDF and making density/CDF/simulation mutually consistent; finite support endpoints are handled accordingly.

## EM organization

- `emstep_ph` / `fit_ph_em`: exact/censored PH EM via Van Loan blocks.
- `emstep_dph` / `fit_dph_em`: discrete PH EM.
- `emstep_bivph`, `emstep_bivdph`: feed-forward bivariate EM.
- `emstep_mdph`: multivariate discrete PH EM.
- `emstep_mph_rc`: right-censored multivariate PH EM.  It uses conditional independence given the shared initial state to avoid the very large temporary tensors in upstream `mPH_EM_UNI.cpp`; the sufficient statistics are algebraically the same.
- `emstep_mphstar`: originating-PH plus marginal-reward decomposition of the upstream MPH* algorithm.

## Not recreated as object APIs

The package's formula/model-frame/S4 machinery, plotting/printing/progress code, `reshape2` reshaping, and `nnet::multinom` mixture-of-experts orchestration are R integration layers and are omitted.  Fixed regression/inhomogeneous likelihood evaluation is available, but the generic `stats::optim` outer loop over transform/regression coefficients is intentionally left to the caller rather than adding another optimizer abstraction.
