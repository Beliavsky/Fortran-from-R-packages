# API

## `bvls_fit`

```fortran
call bvls_fit(a, b, lower, upper, result [, key, istate])
```

Inputs:

- `a(m,n)`: design matrix.
- `b(m)`: response/data vector.
- `lower(n)`, `upper(n)`: component bounds.
- `key`: optional; zero starts from all lower bounds, nonzero uses `istate`.
- `istate(n+1)`: optional warm-start state in the original Stark-Parker form.

The first `istate(n+1)` entries identify bound variables. A positive index is
at its upper bound and a negative index is at its lower bound. Remaining
entries identify active variables. `istate(n+1)` stores the number of bound
variables.

## `bvls_core`

Lower-level interface corresponding closely to the original subroutine:

```fortran
call bvls_core(key, a, b, lower, upper, x, w, istate, loop_a, status)
```

On success, `w(1)` follows the original routine and contains the residual
2-norm. `bvls_fit` exposes the full gradient separately.

## Status constants

- `bvls_success = 0`
- `bvls_max_iterations = 1`
- `bvls_invalid_dimensions = -1`
- `bvls_inconsistent_bounds = -2`
- `bvls_no_free_variables = -3`
- `bvls_invalid_state = -4`
- `bvls_rank_failure = -5`
