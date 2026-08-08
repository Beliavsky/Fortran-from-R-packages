# neighbours-fortran

Modern Fortran/FPM translation of the computational code in Enrico Schumann's
R package `neighbours` 0.1-5.

The library constructs and applies neighbourhood moves for local-search
algorithms. The R package builds closures dynamically; the Fortran translation
stores the same parameters in derived-type configuration objects and applies
explicit procedures.

## Implemented functionality

- numeric zero-sum neighbourhoods (the upstream default `sum = TRUE`)
- numeric one-coordinate moves with no budget (`sum = FALSE`)
- numeric one-coordinate moves with lower/upper total-budget limits
- scalar or per-coordinate lower/upper bounds
- random or fixed step sizes
- active coordinate sets
- incremental `A*x` updating, corresponding to upstream `update = "Ax"`
- logical toggling neighbourhoods
- fixed-cardinality logical neighbourhoods
- bounded-cardinality logical neighbourhoods
- permutation neighbourhoods for real, integer, logical and character vectors
- the upstream `type = "5/10/40"` portfolio neighbourhood
- `next_subset` (Nijenhuis-Wilf NEXKSB ordering)
- the computational part of `compare_vectors`
- upstream internal `random_vector` functionality for numeric/logical vectors

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example numeric_portfolio
fpm run --example logical_fixed_k
```

A compiler-only strict test script is provided in `scripts/` for systems where
FPM is not available.

## NMOF integration

The R package lists NMOF under `Suggests`; it is not a hard dependency. This
translation follows the same design: the core library is standalone.

The user-supplied NMOF Fortran translation is retained under `vendor/nmof`, and
`integration/nmof-demo` is a separate FPM project showing how to pass a
`neighbours` numeric move directly to `nmof_optimization::local_search` via the
exact `real_neighbour` callback signature.

## License

GPL-3.0-only. The original R package is retained under
`original/neighbours-master/`. The supplied NMOF translation retains its own
GPL-3.0-only notices under `vendor/nmof/`.
