# RCEIM-fortran

Modern Fortran translation of the computational portions of the R package **RCEIM 0.3**, with an FPM project layout.

RCEIM implements a stochastic, derivative-free optimizer inspired by the Cross-Entropy Method. The translation keeps the population/elite update equations, smoothing parameters, super-individual carryover, boundary clipping, wait-generation stopping rule, chaos/restart mechanism, minimization/maximization convention, and convergence history.

## Build

```text
fpm build
fpm test
```

Two examples are included:

```text
fpm run --example quadratic_example
fpm run --example rceim_benchmark_example
```

For a strict GNU Fortran build independent of FPM:

```text
scripts/test_gfortran.sh
```

or on Windows:

```text
scripts\test_gfortran.bat
```

## Basic use

```fortran
use rceim, only : dp, rceim_options, rceim_result, ceim_optimize

type(rceim_options) :: opt
type(rceim_result) :: res
real(dp) :: lower(2), upper(2)

lower = [-10.0_dp, -10.0_dp]
upper = [ 10.0_dp,  10.0_dp]
opt%n_total = 500
opt%n_elite = 125
opt%n_super = 1
opt%epsilon = 0.01_dp
opt%max_iter = 60
opt%seed = 2026

call ceim_optimize(objective, lower, upper, res, opt)
```

The objective interface is

```fortran
function objective(x) result(value)
    use rceim, only : dp
    real(dp), intent(in) :: x(:)
    real(dp) :: value
end function objective
```

`res%x` contains the best parameters, `res%value` the objective in its natural sign, and `res%score` the internally minimized score. Thus for maximization, `score = -value`, matching RCEIM's `mfactor` convention while giving Fortran callers the natural objective value as well.

## Scope

Translated:

- `ceimOpt`
- `enforceDomainOnParameters`
- score-based row sorting used by `sortDataFrame`
- `testFunOptimization`
- `testFunOptimization2d`
- convergence histories and final elite population

Not translated:

- `plotEliteDistrib`
- `overPlotErrorPolygon`
- convergence graphics
- interactive `readline()` stepping
- R `parallel::mclapply` orchestration
- generic R data-frame behavior unrelated to optimization

See `TRANSLATION_COVERAGE.md` for exact compatibility notes.

## License

The upstream package declares `License: GPL (>= 2)`. This translation is distributed under **GPL-2.0-or-later** and retains the original package under `original/RCEIM-master/` for provenance. See `LICENSE` and `UPSTREAM_PROVENANCE.md`.
