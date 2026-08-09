# API

## Module `cec2013`

The convenience module re-exports the public API from the data and evaluation
modules.

### Kind

```fortran
integer, parameter :: dp = kind(1.0d0)
```

### `type(cec2013_context)`

Stores one dimension's official shift and rotation data.

```fortran
call ctx%init(n, data_dir, status, message)
call ctx%clear()
logical_value = ctx%initialized()
```

`n` must be one of 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.
`data_dir` must contain `shift_data.txt` and `M_D<n>.txt`.

### `cec2013_evaluate`

```fortran
f = cec2013_evaluate(ctx, problem, x, status)
```

- `problem`: integer 1 through 28.
- `x`: vector of length `ctx%n`.
- `status`: optional integer status.

### `cec2013_evaluate_batch`

```fortran
call cec2013_evaluate_batch(ctx, problem, x, f, status)
```

`x` has shape `(dimension, number_of_points)` and `f` has length
`number_of_points`.

### `cec2013_optimum_value`

```fortran
fopt = cec2013_optimum_value(problem)
```

Returns the published objective offset for problems 1 through 28.

### `cec2013_dimension_supported`

```fortran
if (cec2013_dimension_supported(n)) ...
```

## Status constants

- `CEC2013_OK`
- `CEC2013_BAD_DIMENSION`
- `CEC2013_IO_ERROR`
- `CEC2013_BAD_PROBLEM`
- `CEC2013_BAD_SHAPE`
