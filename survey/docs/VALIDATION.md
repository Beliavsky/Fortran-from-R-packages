# Validation

The translated library is tested with GNU Fortran in strict Fortran 2018 mode:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -ffree-line-length-none
```

The active `minqa` vendor source is compiled with warnings disabled rather than
`-Werror`; translated survey sources and the other active dependencies compile
with warning promotion enabled.

Current test programs cover:

1. Taylor estimators and finite-population corrections.
2. JK1 replicate variance, calibration, and post-stratification.
3. Weighted quantile rules and Gaussian survey GLM.
4. Binomial/logit and Poisson/log survey GLMs against independent references.
5. Wald/Rao-Scott contingency tests and survey IV regression.
6. PPS variance, reliability, and PCA/correlation primitives.
7. Special-function tail probabilities and generic survey MLE.
8. Weighted KM, Cox, log-rank, and parametric survival regression integration.
9. Generic nonlinear survey regression against an independent nonlinear least-squares reference.

A clean strict run should print:

```text
test_chisq_ivreg: PASS
test_estimators: PASS
test_glm_families: PASS
test_multivariate_pps: PASS
test_nls: PASS
test_quantiles_glm: PASS
test_replicates_calibration: PASS
test_special_mle: PASS
test_survival: PASS
```
