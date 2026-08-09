# API

## Module

```fortran
use cec2005benchmark
```

`dp` is the public real kind (`real64`).

## `type(cec2005_context)`

A context caches the data needed for one benchmark function and dimension.

### `init`

```fortran
call ctx%init(function_id, dimension, data_dir, noise_enabled, ierr)
```

- `function_id`: integer 1 through 25.
- `dimension`: one of 2, 10, 30, 50.
- `data_dir`: optional; defaults to `data`.
- `noise_enabled`: optional; defaults to `.true.`.
- `ierr`: optional.  Zero means success; 1 is an invalid function id; 2 is an
  unsupported dimension.  File-I/O errors are returned as compiler I/O status
  values.

### `evaluate`

```fortran
f = ctx%evaluate(x, ierr)
```

Evaluates a single vector.  `x` must have the dimension used to initialize the
context.

### `evaluate_batch`

```fortran
call ctx%evaluate_batch(x, f, ierr)
```

`x(npoints, dimension)` contains one point per row and `f(npoints)` receives
its objective values.

### `set_noise`

```fortran
call ctx%set_noise(.false.)
```

Turns stochastic noise on or off.  For composition functions whose
normalization itself is noisy, normalization constants are recomputed.

## Convenience entry points

```fortran
f = cec2005_eval(function_id, x, data_dir, noise_enabled, ierr)
call cec2005_eval_batch(function_id, x, f, data_dir, noise_enabled, ierr)
```

These construct a temporary context, so a persistent context is preferable in
optimization loops.

## RNG seed

```fortran
call cec2005_seed(12345)
```

Seeds the intrinsic Fortran generator used by the noisy functions.  The noise
distribution follows the C implementation, but the stream does not reproduce
R's RNG bit for bit.
