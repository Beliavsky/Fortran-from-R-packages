# tvGarchKF modern Fortran

A modern Fortran/FPM translation of the computational core of the R package
`tvGarchKF` 0.0.1. The package estimates slowly varying GARCH(1,1)
coefficients through the state-space/Kalman recursion described by the
upstream implementation.

## Implemented functionality

- Polynomial, nonlinear-power, and trigonometric parameter paths
- Built-in `u` and `3*(1-log(u))` trigonometric arguments
- Custom scalar argument functions through a procedure pointer
- The two upstream C Kalman recursions in one array-valued Fortran routine
- Missing observations and appended forecast periods
- The upstream objective used by `tvGarchKalmanLoglike`
- Coefficient fitting and five-decimal compatibility output
- Time-varying GARCH simulation
- Rolling local GARCH(1,1) estimation corresponding to `tvParameter`
- Optional corrected GARCH positivity constraints
- Typed model, filter, fit, simulation, and rolling-estimate results

Plotting, R expression parsing, data frames, `zoo`/time-series classes, and R
printing infrastructure are intentionally omitted.

## Basic use

```fortran
use fgarch_kinds, only : dp
use tvgarchkf

type(tvgarch_spec) :: model
type(tvgarch_simulation_result) :: sim
type(tvgarch_fit_result) :: fit

model = make_tvgarch_spec( &
   make_tv_function([0.04_dp, 0.03_dp], tv_polynomial), &
   make_tv_function([0.08_dp, 0.03_dp], tv_polynomial), &
   make_tv_function([0.78_dp,-0.08_dp], tv_polynomial))

sim = tvgarch_simulate(500, model, seed=1234, &
                       corrected_constraints=.true.)
fit = tvgarch_kalman_fit(sim%returns, model, &
                         corrected_constraints=.true.)
```

For the upstream example argument `3*(1-log(u))`, set
`argument_kind=arg_three_one_minus_log` in a trigonometric specification.
Other expressions are represented by assigning a pure scalar procedure to
`spec%custom_argument` and setting `argument_kind=arg_custom`.

## Build

With FPM:

```text
fpm test
fpm run
```

With GNU Make:

```text
make checked
make optimized
```

The checked build uses Fortran 2018, warnings as errors, bounds/runtime checks,
and backtraces. The optimized build uses `-O3` with the same strict warning
policy.

## Main modules

- `tvgarchkf_types`: typed specifications and results
- `tvgarchkf_functions`: deterministic coefficient paths
- `tvgarchkf_core`: filter, objective, fit, simulation, rolling estimation
- `tvgarchkf`: public umbrella module
- `fgarch_*`: vendored compatible dependency used by `tv_parameter`

See `API_MAP.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.
