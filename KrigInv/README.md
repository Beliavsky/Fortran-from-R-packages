# KrigInv-fortran

Modern Fortran/FPM translation of the computational code in **KrigInv 1.4.2**.

Version **0.2.0** adds the DiceKriging model-estimation layer required to reproduce KrigInv's normal `CovReEstimate=TRUE` sequential behavior.  The project is self-contained: the required subset of the previously translated DiceKriging 1.6.1 computational core is vendored under its GPL-3 license option.

## Included computational functionality

The Fortran port contains:

- `fit_krig_model`, a native DiceKriging-compatible Gaussian-process fitting path;
- covariance/range and variance estimation by MLE, PMLE/SCAD, or LOO where supported by DiceKriging;
- bounded multistart BFGS fitting with analytic MLE/LOO gradients;
- tensor-product and isotropic Gaussian, exponential, Matern 3/2, Matern 5/2, and power-exponential covariance models;
- nonlinear DiceKriging scaling covariances;
- constant, linear, linear-with-interactions, and quadratic trend bases;
- nugget and heteroskedastic observation-noise handling;
- `update_krig_model(..., cov_reestimate=.true.)`, including covariance/variance refitting after new observations;
- the KrigInv-compatible default `CovReEstimate = model@param.estim` in `egi` and `egi_parallel`;
- `init_dice_krig_model` for fixed covariance parameters using DiceKriging covariance conventions;
- the original v0.1.0 `init_krig_model` fixed-parameter path for backward compatibility;
- SK/UK prediction, posterior covariance, and conditional-update algebra;
- Ranjan, Bichon, TMSE, and TSEE pointwise criteria;
- SUR, JN/volume-variance, TIMSE/IMSE, Vorob'ev, conservative-Vorob'ev, and future-volume batch criteria;
- Sobol, Monte Carlo, and KrigInv importance-sampling integration designs;
- discrete and bounded continuous optimization for infill criteria;
- sequential `egi` and batch `egi_parallel` workflows;
- conservative excursion-level calculations used by `vorobCons`/`vorobVol`;
- a native bivariate-normal CDF replacing `pbivnorm`.

Plotting routines (`print_uncertainty*`) and R formula/S4/data-frame/cluster infrastructure are intentionally omitted.

## CovReEstimate behavior

A fitted model records whether DiceKriging parameters were estimated.  For example:

```fortran
use kriginv

type(krig_model) :: model
call fit_krig_model(model,x,y,covariance='matern5_2',trend_order=0)
```

For such a model, the default sequential update is equivalent to KrigInv's:

```r
kmcontrol$CovReEstimate <- model@param.estim
update(..., cov.reestim = kmcontrol$CovReEstimate)
```

Thus:

```fortran
call update_krig_model(model,newx,newy)
```

re-estimates covariance parameters when `model%param_estim` is true.  The behavior can be overridden explicitly:

```fortran
call update_krig_model(model,newx,newy,cov_reestimate=.false.)
```

`egi` and `egi_parallel` accept the same optional `cov_reestimate`, `trend_reestimate`, and `nugget_reestimate` controls.  If `cov_reestimate` is omitted, the default is `model%param_estim`, matching upstream KrigInv.

## Build

```text
fpm build
fpm test
fpm run --example basic_kriginv
fpm run --example cov_reestimate
```

The source uses Fortran 2018 and has no R/C/C++ runtime dependency.

## Licenses

KrigInv is GPL-3; the full license is in `LICENSE`.

The vendored DiceKriging computational modules are from DiceKriging 1.6.1, which upstream distributes under `GPL-2 | GPL-3`.  They are incorporated here under the **GPL-3 option**; provenance and both upstream license texts are retained under `licenses/DiceKriging/`.

The bundled conservative-estimation modules derive from the GPL-3 Fortran `anMC` translation.  Sobol direction-number material retains its randtoolbox/upstream BSD notices under `licenses/`.
