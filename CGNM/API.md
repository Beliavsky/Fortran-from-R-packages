# API

## Main types

`type(cgnm_problem)` stores the model callback, target, initial ranges, bound
masks/values, optional multi-objective parameter targets, and the
`keep_initial_distribution` mask.

`type(cgnm_options)` contains:

- `num_minimizers`
- `num_iterations`
- `initial_lambda`
- `gamma` / `gamma_auto`
- `stay_in_initial_range`
- `seed`
- `max_initial_draws`
- `kmeans_max_iter`

`type(cgnm_result)` returns internal and physical parameters, model values,
initial cluster, residual/lambda histories, target/weight matrices, status, and
iteration count.

## Model callback

```fortran
subroutine model(x, y, ierr)
   use cgnm, only : dp
   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: y(:)
   integer, intent(out) :: ierr
end subroutine model
```

Set `ierr /= 0` when the model cannot be evaluated at a proposed point.

## Problem construction

```fortran
call cgnm_init_problem(prob, model, target, initial_lower, initial_upper, &
   lower_bound=lb, upper_bound=ub, lower_active=has_lb, &
   upper_active=has_ub, mo_weights=mow, mo_values=mov, &
   keep_initial_distribution=keep, ierr=ierr)
```

The logical bound masks are the Fortran equivalent of finite values versus
`NA` in the R interface.

## Main solve

```fortran
call cgnm_fit(prob, opt, result)
call cgnm_fit(prob, opt, result, initial_iterates=theta0)
call cgnm_fit(prob, opt, result, target_matrix=targets)
call cgnm_fit(prob, opt, result, weight_matrix=weights)
call cgnm_fit(prob, opt, result, algorithm_version=1)
```

`initial_iterates` are supplied in physical parameter coordinates. Internally
they are transformed when finite bounds are active.

## Extensions

```fortran
call cgnm_bootstrap(prob, opt, base, nboot, bootstrap_type, boot)
```

`bootstrap_type` values mirror the R package:

1. residual resampling,
2. multinomial observation counts,
3. independent uniform random weights.

```fortran
call cgnm_ebe(prob, opt, base, individual_indices, nrepeat, ebe_weight, ebe)
```

## Postprocessing

```fortran
call top_indices(result, n, idx)
call accepted_indices(result, idx)
call accepted_indices_binary(result, mask)
call best_approximate_minimizers(result, n, theta, idx)
threshold = accepted_max_ssr(result)
call column_quantiles(data, probability, q)
```

The Grubbs-screening path in `accepted_max_ssr` uses a standalone Student-t
quantile approximation rather than R's `qt`.

## Low-level helpers

```fortran
call cgnr_ata_atb(A, b, x)
call cgnr_ata_atb_reg(A, b, lambda, x)
call cgnm_physical_to_internal(prob, theta, z, ierr)
call cgnm_internal_to_physical(prob, z, theta)
call cgnm_evaluate(prob, z, y, ierr)
```
