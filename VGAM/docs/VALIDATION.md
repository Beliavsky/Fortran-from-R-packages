# Validation

Validated on 2026-08-14 with GNU Fortran using:

```text
gfortran -std=f2018 -Werror=implicit-interface -fcheck=all
```

Clean full-tree results:

```text
test_core: PASS
test_models: PASS
test_v02: PASS
test_v03: PASS
test_v04: PASS
test_v05: PASS
test_v06: PASS
test_v07: PASS
test_v08: PASS
test_v09: PASS
```

`test_core` covers reference numerical values and inversion identities for
standard, VGAM-specific, actuarial, and extreme-value distributions and links.

`test_models` covers Gaussian/Poisson IRLS, multinomial and ordinal fitting,
beta regression, negative-binomial regression, zero-inflated Poisson, constrained
multi-response VGLM fitting, and spline-backed GAM fitting.

`test_v02` covers RR-VGLM, zero-truncated/hurdle/ZI count models, finite-support
GAITD transforms, Yeo-Johnson fitting, and Gaussian AR(1).

`test_v03` covers rank-1/rank-2 full symmetric QRR surfaces (including explicit
cross-latent curvature), three-predictor LMS/Yeo-Johnson, and GARMA(p,0).

`test_v04` adds:

- rank-2 DRR recovery with nontrivial `H.A` and `H.C` constraints;
- rank-1 nested RRAR coefficient recovery and recursive forecasting;
- covariate-dependent GAITD inflated-mass regression;
- CQO response-to-latent calibration followed by response-surface reconstruction;
- rank-1 spline CAO with a strong improvement over an intercept-only model.


`test_v05` adds:

- exact VGAM MLM GAITD normalization checks with alteration, inflation,
  additive deflation, and truncation for Poisson and negative-binomial parents;
- synthetic direct-MLM GAITD Poisson regression recovery;
- independence identities for Clayton, Frank, FGM, Gaussian, and Plackett copulas;
- the analytic Gaussian copula quadrant probability at `(u,v,rho)=(0.5,0.5,0.5)`;
- seeded FGM simulation followed by copula maximum-likelihood recovery.


`test_v06` adds:

- Student-t Cauchy-limit density/CDF/quantile identities and direct verification
  of the upstream bivariate Student-t log-density formula;
- an external SciPy grid check of the Student-t CDF, with maximum absolute error
  about `4.3e-16`;
- seeded bivariate Student-t simulation followed by joint df/rho likelihood recovery;
- seeded fixed-df Student-t copula simulation and correlation recovery;
- AMH independence/RNG identities and seeded AMH copula-regression recovery;
- exact GAITD outer-distribution mix allocation checks for alteration, inflation,
  deflation, and truncation, plus a negative-binomial normalization check;
- zero/one-altered beta and zero/one-inflated beta-binomial identities.

`test_v07` adds:

- score outer-product information, constrained information projection, and
  covariance lifting identities;
- censored Poisson/exponential probability identities and seeded censored-normal
  likelihood recovery with simultaneous left and right censoring;
- seeded GAITD outer-mix Poisson simulation followed by parent mean, altered-mass,
  and outer-mean recovery;
- seeded bivariate-normal five-parameter recovery;
- VGAM bivariate-logistic CDF identity plus seeded simulation/likelihood recovery;
- seeded Freund (1961) bivariate exponential four-parameter recovery.


`test_v08` adds:

- positive-normal and positive-geometric identities plus normalization/inversion
  checks for zero-altered Poisson/NB/binomial and zero-deflated geometric laws;
- seeded zero-altered Poisson regression recovery for both parent mean and
  observed-zero probability;
- trivariate-normal density/log-density identity, seeded simulation, and
  nine-parameter mean/scale/correlation likelihood recovery;
- FGM-exponential independence identity and seeded dependence recovery;
- GAITD outer-mix negative-binomial fitting with distinct parent/altered/
  inflated/deflated dispersion predictors on three-point restricted supports.


`test_v09` adds:

- exact Dirichlet density and analytic expected-information identities, plus
  seeded joint Dirichlet-regression shape recovery;
- Tobit endpoint-mass/quantile identities and seeded mean/scale recovery under
  simultaneous lower and upper censoring;
- generalized folded-normal density/quantile identities and seeded likelihood
  recovery for asymmetric fold constants;
- positive negative-binomial normalization/quantile identities and seeded
  parent-mean/dispersion recovery;
- seeded zero/one-altered beta recovery for beta mean, precision, and both
  endpoint masses.

The v0.9 example produces on the validation compiler:

```text
Dirichlet fitted shapes:    2.1424    3.2107    5.3301
Tobit mean coefficients:    0.3260    0.7624
Tobit fitted sd:    0.8237
```

The analytic trivariate-normal log-density was also compared externally with
SciPy over 200 random positive-definite covariance cases; maximum absolute
log-density difference was approximately `1.1e-12`.

The v0.8 example produces on the validation compiler:

```text
Zero-altered Poisson lambda:     2.4562
Observed zero probability:        0.1960
Trivariate normal means:         0.1889   -0.4721    0.7609
Trivariate normal correlations:     0.3861   -0.1636    0.2742
```

The v0.7 example produces on the validation compiler:

```text
Censored-normal fitted mean:    0.4376
Censored-normal fitted sd:      1.1230
Bivariate-normal means:         0.2942  -0.4358
Bivariate-normal sds:           1.1229   0.7497
Bivariate-normal rho:           0.5768
```

The v0.6 example produces on the validation compiler:

```text
GAITD mix P(Y=0):   0.080321
GAITD mix P(Y=4):   0.154381
GAITD mix mean:     2.267181
Bivariate-t fitted df:     7.35251
Bivariate-t fitted rho:    0.48355
```

The v0.5 example demonstrates a direct MLM GAITD distribution and FGM copula
simulation/fitting. On the validation compiler it produces:

```text
GAITD P(Y=0):   0.12000
GAITD P(Y=1):   0.25265
GAITD mean:     2.08265
FGM true parameter:     0.55000
FGM fitted parameter:   0.55321
```

The v0.4 example produces, on the validation compiler:

```text
DRR residual deviance:     0.000000
DRR constrained loadings:     0.8732    1.7464
CAO residual deviance:     0.009596
CAO canonical coefficients:     1.0000    0.5489
RRAR fitted lag-1 matrix:
    0.5573    0.1347
    0.2570    0.0621
```

All `.f90` source, test, and example lines are checked to be no longer than 132
characters. The standard FPM source/test/example layout does not require
compiler-specific line-length flags.
