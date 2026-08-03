# API

All public entities are available through:

```fortran
use segmented
```

All real calculations use `dp = kind(1.0d0)`.

## Controls and results

### `type(segmented_control)`

Important fields:

- `max_iter`: outer breakpoint iterations.
- `glm_max_iter`: IRLS iterations for each GLM fit.
- `grid_points`: maximum threshold intervals inspected per step breakpoint.
- `max_line_search`: safeguarded continuous-breakpoint line-search iterations.
- `tolerance`: regression/IRLS tolerance.
- `breakpoint_tolerance`: relative breakpoint convergence tolerance.
- `lower_quantile`, `upper_quantile`: admissible breakpoint range.
- `verbose`: iteration output.

### `type(segmented_result)`

Contains `coefficients`, `covariance`, `breakpoints`, `breakpoint_se`, `fitted`,
`residuals`, `weights`, objective and likelihood statistics, AIC/BIC, iteration
history, status, and convergence flag.

### `type(segmented_lme_result)`

Contains an `nlme` `lme_result` in `fit`, breakpoint estimates and standard
errors, history, status, and convergence flag.

## Fitting routines

```fortran
call fit_segmented_lm(y, x, z, psi0, result, weights, offset, control)
call fit_stepmented_lm(y, x, z, psi0, result, weights, offset, control)
call fit_segmented_glm(y, x, z, psi0, family, result, weights, offset, control)
call fit_stepmented_glm(y, x, z, psi0, family, result, weights, offset, control)
call fit_segmented_lme(y, x, z, psi0, random_design, group, result, options, control)
```

Aliases are provided as `segmented_lm`, `stepmented_lm`, `segmented_glm`,
`stepmented_glm`, and `segmented_lme`.

`family` is one of:

- `FAMILY_GAUSSIAN`
- `FAMILY_BINOMIAL`
- `FAMILY_POISSON`

`x(n,p)` is the ordinary fixed-effect matrix. `z(n,m)` holds the covariate for
each estimated breakpoint. To estimate several breakpoints on one covariate,
repeat that covariate in several columns of `z`.

For continuous fits the returned coefficient order is:

```text
base X coefficients, hinge/slope-change coefficients
```

For step fits it is:

```text
base X coefficients, post-break level changes
```

## Convenience interfaces

```fortran
call segmented_numeric(y, result, x_values, psi0, n_break, weights, control)
call stepmented_numeric(y, result, x_values, psi0, n_break, weights, control)
call segreg(y, x, z, psi0, family, kind, result, weights, offset, control)
call stepreg(y, x, z, psi0, family, result, weights, offset, control)
```

## Model matrices and prediction

```fortran
call segmented_model_matrix(x, z, breakpoints, kind, design, status)
call predict_segmented(fit, xnew, znew, prediction, status, linear_predictor)
call hinge_matrix(z, breakpoints, u, v)
call step_matrix(z, breakpoints, h)
```

`v` is the negative indicator derivative used by the Muggeo linearization.

## Interpretation and inference

```fortran
call segment_slopes(fit, base_coefficient, slopes, standard_errors, status)
call segment_intercepts(fit, intercept_coefficient, base_slope_coefficient, &
    intercepts, status)
call breakpoint_confint(fit, confidence, interval, status)
value = aapc(fit, x_min, x_max, base_slope_coefficient, exponentiate, status)
y = broken_line_values(x, intercept, base_slope, changes, breakpoints)
```

The slope/intercept helpers assume that the breakpoint columns refer to one
ordered segmented covariate.

## Tests, power, and selection

```fortran
call davies_test(y, x, z, grid_points, result, weights, offset)
call pscore_test(y, x, z, breakpoint, result, weights, offset)
power = pwr_seg(z, x, breakpoint, slope_change, sigma, alpha, two_sided, status)
call select_breakpoints_bic(y, x, z, max_breaks, family, best_fit, bic_values, &
    weights, offset, control)
```

`davies_test` uses a conservative grid/Bonferroni calibration. `pwr_seg` is a
fixed-breakpoint normal-score approximation.
