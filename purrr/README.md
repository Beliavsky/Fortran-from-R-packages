# purrr

This package is a restricted modern free-form Fortran translation of
**purrr 1.2.2**, *Functional Programming Tools*. It provides typed functional
operations over intrinsic real and integer arrays using explicit callback
interfaces checked by the compiler.

Implemented operations include typed and generic `map`, `map2`, `modify`,
`reduce`, `accumulate`, `keep`, `discard`, `some`, `every`, `none`, `detect`,
`detect_index`, `head_while`, `tail_while`, and side-effecting `walk`.

```fortran
values = map([1, 2, 3], square)
total = reduce(values, add)
```

Pure transformations require pure scalar callbacks. `walk` deliberately
accepts an impure callback for output or other side effects.

## Build and test

```text
fpm build
fpm test
fpm run --example map_reduce
```

## Scope

The translation uses homogeneous typed arrays rather than R's heterogeneous
lists. It does not reproduce formulas, dynamic dots, environments, arbitrary
object dispatch, asynchronous execution, or condition capture. See
`API_COVERAGE.md` for exact mappings.

## License

MIT, matching upstream purrr. See `LICENSE`, `NOTICE.md`, and `upstream/`.
