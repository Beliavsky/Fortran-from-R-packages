# Porting notes

## Numerical dependencies

The translation intentionally uses the supplied dependency ports.

`survival_f90` supplies the Cox Newton solver and right-censored concordance
calculation used by the Wald/CV layers. Only the modules needed by compound.Cox
are packaged under `vendor/survival-fortran`.

`numDeriv-fortran` supplies Richardson/Bates-Watts Hessians for the variance
calculations in `compound_reg` and `depend_cox_reg`.

## Optimizer replacement

Upstream uses R's `nlm` for its package-specific objectives. The Fortran port
uses a deterministic BFGS optimizer with central numerical gradients. The final
Hessian for standard errors is then computed by `numDeriv-fortran`.

This keeps the same likelihood/objective functions without requiring an R
optimizer runtime. Iteration paths are not expected to be bit-for-bit identical
to `nlm`.

## Cox likelihood conventions

The custom likelihoods in `compound.reg` are translated exactly as risk-set
sums (Breslow-style handling at tied times). `uni_wald` uses the supplied
survival translation's Efron Cox fit, matching the default `survival::coxph`
behavior used by upstream.

## Dependent censoring

`depend_cox_reg` preserves the upstream Clayton-copula joint likelihood. Baseline
survival and censoring cumulative hazard increments are parameterized on the log
scale exactly as in the R source. The implementation evaluates
`log(exp(a)+exp(b)-1)` in a stable log-sum form.

## Compound regression CV

The shrinkage criterion is the upstream convex combination of the multivariate
Cox partial log likelihood and the sum of univariate partial log likelihoods.
The grid search and Verweij-style cross-validated likelihood follow the R code.
The variance correction uses the upstream `h_dot`, `E_z`, `A_hat`, and
`Sigma_hat` formulas.

## Copula-graphic estimators

The Clayton, Frank and Gumbel recursions are direct translations. Frank's
Kendall tau requires the Debye integral; the Fortran version evaluates it with a
high-order composite Simpson rule and a local series for `x/(exp(x)-1)` near
zero.

## Factorial survival analysis

The copula-graphic pairwise effect estimates, delete-one jackknife covariance,
contrast projection, simulated weighted-chi-square critical values, and
Satterthwaite analytical approximation are translated. `MASS::ginv` is replaced
by a symmetric eigenvalue pseudoinverse.

## R-only features

Plotting and R object/data-set handling are excluded. The `randomize=TRUE`
convenience option from `compound.reg` / `uni.selection` is represented by
shuffling input rows before calling the Fortran procedure; default sequential
folds match upstream `randomize=FALSE`.
