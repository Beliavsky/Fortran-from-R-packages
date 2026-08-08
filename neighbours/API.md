# API

All public functionality is provided by module `neighbours`, except the kind
parameters in `neighbours_kinds`.

## Kinds and RNG

- `dp = kind(1.0d0)`
- `i8`
- `type(rng_state)`
- `rng_seed`, `rng_uniform`, `rng_integer`, `rng_shuffle`

The deterministic RNG is adapted from the supplied NMOF Fortran translation so
examples/tests do not depend on processor-global `random_number` state.

## Numeric neighbourhoods

`type(numeric_neighbour_config)` stores bounds, active coordinates, step size,
budget mode, and optional `A` for incremental matrix-product updates.

Initialize with:

```fortran
call init_numeric_neighbour(config, n, stepsize, lower=[0.0_dp], &
     upper=[1.0_dp], status=status)
```

Budget constants:

- `budget_zero_sum`: two-coordinate transfer; total sum unchanged
- `budget_none`: one coordinate moves up or down
- `budget_range`: one coordinate moves subject to total lower/upper budget

Apply with:

```fortran
call numeric_neighbour(config, x, xn, rng, status=status)
```

For upstream `update="Ax"` behavior:

```fortran
call init_numeric_neighbour(config, n, step, a=A, status=status)
Ax = matmul(A, x)
call numeric_neighbour(config, x, xn, rng, ax=Ax, status=status)
```

`Ax` is updated by `A*(xn-x)` rather than recomputed.

## Logical neighbourhoods

Use `type(logical_neighbour_config)`, `init_logical_neighbour`, and
`logical_neighbour`.

With no `kmin/kmax`, `stepsize` distinct entries are toggled. With
`kmin==kmax`, equal numbers of TRUE and FALSE entries are exchanged so the
cardinality stays constant. With `kmin<kmax`, the upstream bounded-cardinality
logic is used.

## Permutations

Generic `permute_neighbour` supports arrays of:

- `real(dp)`
- `integer`
- `logical`
- fixed-length `character`

A missing/one step size corresponds to the upstream default two-position swap.

## Other routines

- `portfolio_5_10_40_neighbour`
- `random_numeric_vector`, `random_numeric_vectors`
- `random_logical_vector`, `random_logical_vectors`
- `next_subset(a,n,k,has_next,status)`
- `compare_logical_vectors(vectors,distances,difference_mask,status)`

Status values are `neighbours_ok`, `neighbours_invalid_input`, and
`neighbours_no_candidate`.
