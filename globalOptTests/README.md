# globalOptTests-fortran

Modern Fortran translation of the computational code in the R package
`globalOptTests` 1.1.

The library provides all 50 benchmark objective functions from the upstream
`src/objFun.c`, plus array-level equivalents of the R helpers:

- `go_test(x, fn_name)`
- `get_default_bounds(fn_name, lower, upper)`
- `get_problem_dimension(fn_name)`
- `get_global_opt(fn_name)`

The individual benchmark functions are also public Fortran procedures.

## Build

```text
fpm build
fpm test
```

No external numerical libraries are required.

## Example

```fortran
use global_opt_tests, only : dp, go_test
real(dp) :: x(2), f

x = [-3.14159265359_dp, 12.275_dp]
f = go_test(x, 'Branin')
```

## Translation policy

The numerical formulas and R metadata are preserved as closely as possible.
Two cases of undefined C behavior are made safe without changing the intended
mathematics: Hartman3/Hartman6 use their four defined table rows, and Gulf
skips the `j=0` evaluation of `log(0)` whose limiting contribution is zero.
Other upstream oddities are preserved and documented in
`TRANSLATION_COVERAGE.md`.

## License

The upstream package declares `GPL (>= 3)`.  The translated code is therefore
provided under GPL-3.0-or-later.  See `LICENSE` and `UPSTREAM_PROVENANCE.md`.
