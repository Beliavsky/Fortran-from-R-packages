# Algorithm and translation notes

## v0.2.0: DiceKriging-backed model fitting

KrigInv itself delegates Gaussian-process fitting to `DiceKriging::km`.  Version 0.1.0 translated KrigInv's criteria and prediction/update algebra but kept covariance hyperparameters fixed.  Version 0.2.0 vendors the numerical core of the separate DiceKriging-fortran translation so that this dependency boundary is now covered inside the FPM project.

`krig_model` has two modes:

1. **Legacy fixed mode**, created by `init_krig_model`, preserving the v0.1.0 API and covariance conventions.
2. **DiceKriging-backed mode**, created by `fit_krig_model` or `init_dice_krig_model`.  Prediction and covariance evaluation use the DiceKriging covariance parameterization, and the model retains the full fitting state needed for re-estimation.

The DiceKriging-backed mode supports Gaussian, exponential, Matern 3/2, Matern 5/2, and power-exponential covariance families, tensor-product or isotropic ranges, scaling covariances, nugget/noise, standard DiceKriging trend bases, MLE/PMLE/LOO fitting, and analytic likelihood/LOO gradients.

## CovReEstimate semantics

Upstream KrigInv sets:

```r
if (is.null(kmcontrol$CovReEstimate))
    kmcontrol$CovReEstimate <- model@param.estim
```

and passes that value to `DiceKriging::update`.  Version 0.2.0 reproduces this behavior:

- a fitted model records `param_estim` from the DiceKriging backend;
- `cov_reestimate_default` is initialized from `param_estim`;
- `update_krig_model` uses that default unless `cov_reestimate` is supplied explicitly;
- `egi` and `egi_parallel` use the same default and forward covariance, trend, and nugget re-estimation controls to each model update.

As in DiceKriging, trend re-estimation defaults to true and nugget re-estimation defaults to false.  When covariance re-estimation is disabled, the current covariance and process variance remain fixed while the trend may still be recomputed.

A zero `new_noise` supplied by EGI does not unnecessarily convert an otherwise noise-free model into the heteroskedastic-noise fitting case.  A genuinely positive new noise variance does activate the noisy model and is included in the refit.

## Kriging algebra

For DiceKriging-backed models, SK/UK prediction is delegated to the translated DiceKriging kernel, including its covariance conventions and universal-kriging correction.  KrigInv's public prediction structure is populated from that result.  Posterior cross-covariance uses the same backend covariance calculation, so criteria and quick-update calculations remain consistent with the fitted model.

For legacy fixed models, the v0.1.0 native algebra remains available unchanged.

## Optimizer replacement

DiceKriging's R `optim` path and KrigInv's external `rgenoud` dependency are not runtime dependencies of this port.  DiceKriging parameter fitting uses the translated bounded multistart BFGS implementation; KrigInv criterion optimization continues to use bounded differential evolution or the translated discrete-search path.

## Conservative excursion sets

The `anMC::conservativeEstimate` dependency used by `vorobCons` and `vorobVol` remains implemented by the bundled GPL-3 modern-Fortran `anMC` modules.

## Intentional corrections retained from v0.1.0

1. **Prospective observation noise.** Requested new-observation variance is added to the prospective covariance diagonal in batch conditioning.
2. **Vorob'ev threshold boundary.** Order-statistic interpolation is clamped to valid first/last indices.
3. **JN importance weights with zero mass.** Uniform proposal mass is used when excursion-probability mass is zero.
4. **Singular numerical designs.** The legacy fixed model gets one small covariance-scaled Cholesky jitter retry.

## Unused upstream import

KrigInv declares `mvtnorm::pmvnorm` in `NAMESPACE`, but the supplied R source does not call it.  No `mvtnorm` runtime dependency is required.

## Deliberately not translated

The `print_uncertainty*` functions are plotting code.  R formula parsing/model frames, S4/S3 display machinery, data-frame/name validation, cluster management, and console progress output are interface infrastructure and are omitted.
