# Porting notes

## Data orientation

Vectors are rank-one arrays. Forecast matrices use observations by models. The returned matrices have the following columns:

1. actual;
2. sGARCH;
3. GJR-GARCH;
4. EWMA/iGARCH;
5. MEM;
6. ensemble.

The accuracy matrix has five rows in the same model order, excluding `actual`. Columns 1-9 are training metrics and columns 10-18 are test metrics.

## Split behavior

The training length is `int(n * split_ratio)`, bounded to leave at least one test observation. This corresponds to the practical truncation used by the R indexing expression.

## GARCH output

The R package uses `model@fit$fitted.values` and `forecast@forecast$seriesFor`. Those are mean forecasts, even though the package describes a volatility ensemble. `PWEV_GARCH_MEAN` is therefore the compatibility default. `PWEV_GARCH_SIGMA` is supplied as an explicit corrected/useful alternative.

## MEM out-of-sample behavior

The original call fits with an `out_of_sample` count and then evaluates `cond_vol` on the held-out dependent-variable vector. Thus later held-out predictions use preceding held-out observations. `PWEV_MEM_UPSTREAM_OOS` preserves this behavior. `PWEV_MEM_RECURSIVE_OOS` substitutes prior predictions after the training endpoint.

## Particle swarm optimization

The PSO uses the upstream bounds `[0,1]`, inertia `0.729`, individual and group constants `1.49445`, maximum velocity `2`, training sample size as the default swarm size, and 1000 iterations. The least-squares objective is evaluated through precomputed cross-products, which is algebraically equivalent and much faster than recomputing residuals for every coordinate update.

## Failure handling

A failed base model is replaced by the training sample mean unless `fail_on_base_model_error` is enabled. The final status remains `PWEV_MODEL_FAILURE` so the fallback cannot be mistaken for a successful fit.

## R-specific omissions

R time-series classes, `xts`/`zoo` dates, column names, row names, list dispatch, console progress bars, and printing methods are not reproduced.
