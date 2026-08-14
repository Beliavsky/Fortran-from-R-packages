# mlrMBO-fortran

Modern Fortran/FPM translation of the computational core of R package
`mlrMBO` 1.1.6.

The library implements model-based/Bayesian optimization without an R
runtime. A target function is supplied as an ordinary Fortran procedure
callback. The native surrogate is the modern Fortran DiceKriging Gaussian
process included as an FPM path dependency.

## Implemented numerical functionality

- initial Latin-hypercube and random designs;
- continuous, integer and categorical parameters, plus simple categorical
  equality conditions;
- optimization paths, continuation and termination by iteration/evaluation
  budget or target value;
- mean-response, standard-error, EI, CB, adaptive-CB, AEI and EQI infill
  criteria with mlrMBO's internal minimization sign convention;
- focus-search infill optimization with independent restarts and progressive
  region shrinking;
- filtering/replacement of proposed points that are too close to existing
  design points;
- constant-liar and parallel-CB multipoint proposals;
- direct indicator-based multi-objective optimization using SMS or additive
  epsilon indicators;
- exact non-dominated sorting, reference points, arbitrary-dimensional
  dominated hypervolume and hypervolume contributions;
- ParEGO normalization, discrete grid weights and augmented Tchebycheff
  scalarization;
- an MSPOT-style native candidate-pool implementation using surrogate
  confidence-bound fronts and removal hypervolume contributions;
- GP covariance re-estimation after model updates through DiceKriging;
- log/square-root response transformation helpers.

All Fortran sources are free source form. The project manifest sets
`implicit-typing = false` and `implicit-external = false`.

## Minimal example

```fortran
use mlrmbo

type(mbo_space) :: space
type(mbo_control) :: control
type(mbo_result) :: result

call init_space(space, [mbo_real,mbo_real], [-5.0_dp,0.0_dp], &
  [10.0_dp,15.0_dp])
call init_control(control, 1, [.true.])
control%infill_criterion = crit_ei
control%max_iter = 10
call mbo(space, objective, control, result)
```

See `example/branin_mbo.f90` for a complete program.

## Deliberate boundaries

`mlrMBO` is partly a framework around algorithms supplied by other R
packages. This port does not claim those adapters as mlrMBO algorithms.
External `mlr` learner adapters, CMA-ES, rgenoud, `emoa` NSGA-II, graphics,
parallelMap, data.table/ParamHelpers object plumbing and R error-handling
wrappers are omitted. The Fortran library uses DiceKriging as its native
surrogate and focus search as its primary infill optimizer.

The upstream MOI-MBO batch method is also deferred because it is specifically
implemented using `emoa` SBX/polynomial-mutation evolutionary operators.
Selecting the reserved `batch_moimbo` value therefore raises an explicit
error instead of silently substituting another algorithm.

See `API_MAPPING.md`, `ALGORITHM_NOTES.md` and `VALIDATION.md` for details.
