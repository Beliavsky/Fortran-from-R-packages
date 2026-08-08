# Fortran API

All public symbols are re-exported by module `calibrar`.

## Callback interfaces

```fortran
function scalar_objective(x) result(f)
  real(dp), intent(in) :: x(:)
  real(dp) :: f
end function

subroutine gradient_callback(x, g)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: g(:)
end subroutine

subroutine vector_objective(x, f)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: f(:)
end subroutine
```

Objective invocation is performed at module-level explicit-interface call sites; the library does not rely on nested host-associated procedure callbacks.

## Optimization

- `optim2(par, fn, result, method, lower, upper, options, gr, active)`
- `optimh(...)` -- same standalone interface; retained for calibrar naming compatibility.
- `calibrate(par, fn, phases, result, method, lower, upper, options, gr, replicates)`
- `calibrate_multi(par, fn, nvar, phases, result, lower, upper, options, weights, replicates)`
- `ahres(par, fn, nvar, result, lower, upper, options, weights, active)`
- `ahres_scalar(...)`

`optim_options` contains iteration/evaluation limits, tolerances, stochastic replicate count, objective scaling/maximization controls, RNG seed, and numerical-gradient method.

`ahr_options` contains AHR-ES population size, selection fraction, smoothing coefficient, step size, stopping-rule settings, stochastic replicate count, and RNG seed.

## Numerical gradients

- `numerical_gradient`
- `gradient_forward`
- `gradient_backward`
- `gradient_central`
- `gradient_richardson`
- `gradient_options`

Richardson extrapolation follows the R source's factor-4 extrapolation sequence.

## Fitness/objective utilities

- `fitness_norm2`
- `fitness_lnorm2`
- `fitness_lnorm3`
- `fitness_lnorm4`
- `fitness_lnorm4b`
- `fitness_rangeq`
- `fitness_pois`
- `fitness_normp`
- `fitness_penalty`
- `fitness_multinom`
- `weighted_sum_fitness`
- `calibration_term`
- `calibration_objective_value`

`calibration_term` represents one calibrated output in a flat observed/simulated vector using `first:last`, a fitness type, weight, and active flag. This replaces the R data-frame/list closure representation while preserving the numerical aggregation step.

## Random/statistical utilities

- `set_seed`
- `rand_uniform`
- `rand_normal`
- `rand_exponential`
- `rtnorm_sample`
- `rtnorm_matrix`
- `dmvnorm_pdf`
- `gaussian_kernel_2d`

## Splines

- `spline_result`
- `spline_par`
- `cubic_spline_eval`

The implementation uses cubic not-a-knot interpolation, corresponding to the endpoint behavior intended by R's default `splinefun`/FMM spline for ordinary cases.

## Stopping rules

- `smooth_stop2`
- `smooth_stop3`
- `smooth_stop4`
- `n_stop`
