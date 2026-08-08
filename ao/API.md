# Fortran API

## Main types

`type(ao_options)` controls the AO process. Important fields are:

- `partition`: `AO_PARTITION_SEQUENTIAL`, `AO_PARTITION_RANDOM`,
  `AO_PARTITION_NONE`, or `AO_PARTITION_CUSTOM`;
- `new_block_probability`, `minimum_block_number`;
- `minimize`;
- `iteration_limit`, `seconds_limit`;
- `tolerance_value`, `tolerance_parameter`, `tolerance_history`;
- `base_optimizer`: `AO_BASE_BFGS`, `AO_BASE_NELDER_MEAD`, or
  `AO_BASE_NEWTON`;
- `base_max_iterations` (default 10, matching the R package's default block
  iteration limit);
- `custom_partition(:)` for explicit blocks.

`type(ao_result)` contains `estimate`, `value`, `seconds`, `iterations`,
`converged`, `stopping_reason`, and `details`.

`details` is an `ao_history` with arrays:

```fortran
iteration(:)
value(:)
parameter(:, :)
active(:, :)
seconds(:)
```

Only entries `1:details%n` are populated.

## Objective callbacks

```fortran
function objective(x) result(value)
  real(dp), intent(in) :: x(:)
  real(dp) :: value
end function

subroutine gradient(x, g)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: g(:)
end subroutine

subroutine hessian(x, h)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: h(:,:)
end subroutine
```

Run without derivatives:

```fortran
call ao_optimize(objective, initial, result, options, lower, upper)
```

With an analytic gradient:

```fortran
call ao_optimize(objective, initial, result, options, lower, upper, &
                 gradient=gradient)
```

With gradient and Hessian:

```fortran
options%base_optimizer = AO_BASE_NEWTON
call ao_optimize(objective, initial, result, options, gradient=gradient, &
                 hessian=hessian)
```

The custom parameter-distance callback has signature
`real(dp) function norm(x,y)` and is passed as `parameter_norm=`.

## Custom partitions

```fortran
options%partition = AO_PARTITION_CUSTOM
allocate(options%custom_partition(2))
allocate(options%custom_partition(1)%index(2))
options%custom_partition(1)%index = [1,2]
allocate(options%custom_partition(2)%index(2))
options%custom_partition(2)%index = [2,3]   ! overlaps are allowed
```

## Multiple processes

`ao_optimize_multiple` forms the Cartesian product of columns of `initials`,
`partition_modes`, and `base_methods`, runs them serially, and returns every
process plus the best one. The R package can parallelize these independent
processes with `future.apply`; parallel scheduling is intentionally not part of
the numerical Fortran core.

## Helpers

- `ao_seed(seed)`
- `generate_random_partition(npar,p,minimum_blocks,blocks)`
- `split_estimate(estimate,npar,parts[,ok])`
- `euclidean_parameter_norm(x,y)`
