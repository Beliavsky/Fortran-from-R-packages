# API

## Main types

### `piqp_settings_type`

Fortran derived-type counterpart of `piqp_settings()`. It includes the upstream
fields `rho_init`, `delta_init`, primal/dual tolerances, regularization limits,
iteration limits, interior step `tau`, preconditioner settings, iterative
refinement settings, and verbosity/timing flags.

### `piqp_result_type`

Contains

* `x`, equality dual `y`;
* inequality lower/upper duals `z_l`, `z_u`;
* box lower/upper duals `z_bl`, `z_bu`;
* corresponding slacks `s_l`, `s_u`, `s_bl`, `s_bu`;
* `info` (`piqp_info_type`) with status, iteration count, residuals, objective
  values, duality gap, proximal parameters and step lengths.

### `piqp_model_type`

Persistent dense model with type-bound procedures:

* `setup`
* `solve`
* `update`
* `get_dims`

Problem dimensions are fixed after `setup`.

## One-shot solve

```fortran
call solve_piqp(pmat, c, result, amat, b, gmat, h_l, h_u, x_l, x_u, settings)
```

All arguments after `result` are optional. `P` may also be omitted by keyword
for a linear objective:

```fortran
call solve_piqp(c=c, result=result, x_l=xl, x_u=xu)
```

The generic interface has a second overload accepting Matrix-fortran
`csc_matrix` objects for `P`, `A`, and `G`.

## Status constants

* `PIQP_SOLVED = 1`
* `PIQP_MAX_ITER_REACHED = -1`
* `PIQP_PRIMAL_INFEASIBLE = -2`
* `PIQP_DUAL_INFEASIBLE = -3`
* `PIQP_NUMERICS = -8`
* `PIQP_UNSOLVED = -9`
* `PIQP_INVALID_SETTINGS = -10`

`status_to_string()` returns the compact PIQP status text;
`status_description()` follows the R package's longer descriptions.
