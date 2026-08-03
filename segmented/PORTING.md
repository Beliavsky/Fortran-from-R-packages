# Porting notes

## R objects and formulas

R formulas, model frames, special `seg()` terms, S3 methods, and mutable
`lm`/`glm`/`lme` objects are replaced with explicit arrays and derived types.
This avoids hidden evaluation rules and gives compile-time dimension checking.

## Continuous segmented models

For breakpoint `psi`, the implementation forms

```text
U = max(z - psi, 0)
V = -I(z > psi)
```

and fits the expanded linearized model `[X,U,V]`. If `gamma` and `delta` are the
coefficients of `U` and `V`, the proposed breakpoint is

```text
psi_new = psi + delta / gamma
```

The update is quantile-clamped, ordered when columns of `z` are identical, and
accepted through objective-decreasing step halving. The final reported model is
`[X,U]`.

## Stepmented models

A step effect is discontinuous in the breakpoint, so derivative linearization
is inappropriate. The Fortran port performs deterministic coordinate searches
at midpoints between observed sorted covariate values. This preserves the exact
least-squares/GLM threshold objective and avoids smoothing the indicator.

## Generalized linear models

Binomial-logit and Poisson-log fits use self-contained IRLS. Binomial responses
may be zero/one values or proportions in `[0,1]`. The package's broader ability
to operate on arbitrary R `family` objects is not reproducible without an R
runtime.

## Mixed effects and `nlme`

The project includes the computational modules from the previous `nlme-fortran`
translation. `fit_segmented_lme` alternates the same hinge linearization with
Gaussian LME fits and uses `nlme` likelihoods, covariance structures, random
covariance parameterizations, and BLUPs.

The R package can put change-point parameters themselves into sophisticated
random-effect formulas. This port currently estimates common fixed breakpoints;
users may include random intercepts and slopes through `random_design`, but a
fully nonlinear random breakpoint is not represented.

## Inference

Breakpoint standard errors use the linearized `V` coefficient divided by the
hinge coefficient. Confidence intervals are normal approximations. The Davies
routine is a grid supremum test with conservative Bonferroni calibration rather
than the package's complete Davies bound implementation. BIC selection is
provided for zero through `max_breaks` continuous breakpoints.

## Features not reproduced exactly

- Formula/S3 methods and plotting.
- Bootstrap restart and random restart history.
- Constrained slope matrices and by-variable formula expansion.
- Segmented ARIMA fitting.
- Full profile/score/bootstrap confidence-interval variants.
- Full random-change-point mixed-effects parameterization.
- Selection procedures that mutate R model formulas.

The original R sources are retained under `original/R` for reference.
