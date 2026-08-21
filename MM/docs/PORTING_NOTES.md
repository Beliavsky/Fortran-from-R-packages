# Porting notes

## Parameter representation

The R S4 class `paras` stores a single numeric vector containing the first
`k-1` multinomial probabilities followed by the upper-triangular interaction
parameters.  `paras_type%vals` deliberately keeps this representation so the
mapping is transparent:

1. entries `1:k-1` are `p(1:k-1)`;
2. entries `k:` are `theta(i,j)` in R's upper-triangle column order:
   `(1,2), (1,3), (2,3), ...`.

The final probability is reconstructed as `1-sum(p(1:k-1))` just as upstream.

## Normalizing constant

Upstream computes `NormC` by enumerating all weak compositions and summing
`MM_single` values.  The Fortran version uses the same support from the
supplied `partitions` translation but sums in log space.  This preserves the
model while improving numerical range.

## Lindsey's Poisson device

The R code fits

```text
d ~ -1 + offset(-sum(lfactorial(y))) + (.)^2
```

with a Poisson GLM.  The Fortran design matrix contains the same columns:
all `k` main count variables followed by all pair products `y_i*y_j` in upper
triangle order.  `poisson_glm_fit` maximizes the Poisson log likelihood with
Newton/IRLS steps and backtracking.  Exponentiated main effects are normalized
to produce `p`; exponentiated pair coefficients produce `theta`.

## Likelihood optimizer

Upstream optimizes `log(getVals(paras))`, exponentiates inside the objective,
and reconstructs the final probability as one minus the first `k-1` entries.
The Fortran optimizer retains that transformed parameterization and invalid
simplex/line-search points receive a large objective penalty.

The R package offers `nlm` and Nelder-Mead.  There is no R `nlm` runtime in
this port, so `method="nlm"` (and the default) uses a native finite-difference
BFGS optimizer with Armijo backtracking.  `method="Nelder"` uses a direct
Nelder-Mead implementation.

## R object infrastructure

S3/S4 methods, formula parsing, printing, and Oarray conversions are language
infrastructure rather than numerical algorithms.  They are not reproduced.
`gunter_type` and `gunter_mb_type` expose the computational support tables and
frequency vectors directly.

## Multiplicative-binomial code

`Lindsey_MB` is translated for the same bivariate case as upstream.  Its six
Poisson-regression design columns are intercept, `x1`, `x2`,
`x1*(m1-x1)`, `x2*(m2-x2)`, and `x1*x2`, with the two binomial log
coefficients as the offset.
