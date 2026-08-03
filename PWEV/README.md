# PWEV modern Fortran translation

This project translates the computational workflow of the R package `PWEV` 0.1.0 to modern Fortran and packages it for FPM.

The main entry point is:

```fortran
use pwev

type(pwev_result) :: result
integer :: status

call pwev_fit(data, 0.90_dp, result, status)
```

`result` contains:

- `train_fitted`: actual training data, four base-model fitted series, and the ensemble;
- `test_pred`: actual held-out data, four base-model predictions, and the ensemble;
- `accuracy`: five rows by 18 columns, containing nine training and nine test metrics;
- `weights`: the four PSO ensemble weights;
- model and optimizer status information.

## Base models

The workflow fits:

1. standard GARCH(1,1);
2. GJR-GARCH(1,1);
3. integrated GARCH(1,1), used by the upstream package as EWMA;
4. a no-skew multiplicative error model.

The attached `rugarch` and `rumidas` Fortran translations are vendored as FPM path dependencies.

## Ensemble optimizer

The original package delegates to `WeightedEnsemble`, whose default optimizer is particle swarm optimization. This port implements the same box-constrained objective directly:

```text
minimize sum((actual - forecasts * weights)^2)
subject to 0 <= weight(j) <= 1
```

The weights are not forced to sum to one, matching the upstream R code. The default PSO constants and update order follow the `metaheuristicOpt` implementation used upstream.

## Important compatibility choices

By default, the GARCH columns reproduce the upstream use of `fitted.values` and `seriesFor`, which are conditional-mean outputs. Set `control%garch_output = PWEV_GARCH_SIGMA` to ensemble conditional standard deviations instead.

The upstream MEM out-of-sample calculation evaluates the fitted recursion on the held-out observations themselves. This behavior is retained by default with `PWEV_MEM_UPSTREAM_OOS`. Set `PWEV_MEM_RECURSIVE_OOS` for a genuine recursive forecast that does not use held-out observations as lags.

The unused `rv` expression in the original R function is intentionally omitted.

## Build

```sh
fpm build
fpm test
fpm run
```

GNU Fortran scripts are also supplied:

```sh
scripts/test_gfortran.sh
scripts/test_gfortran_optimized.sh
```

Windows:

```bat
scripts\test_gfortran.bat
scripts\test_gfortran_optimized.bat
```

## License

The combined package is distributed under GPL-3.0-only. See `LICENSE`, `NOTICE.md`, and `THIRD_PARTY_LICENSES.md`.
