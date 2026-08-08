# OOR-fortran

Modern Fortran/FPM translation of the computational code in the R package
**OOR 0.1.4 (Optimistic Optimization in R)**.

The package implements optimistic black-box global optimization for functions
whose gradients may be unavailable or nonexistent.

## Translated algorithms

- `poo`: Parallel Optimistic Optimization (POO), using the package's shared
  HOO-style binary tree.
- `stosoo`: stochastic StoSOO and deterministic SOO in arbitrary dimension.
- Benchmark functions: `guirland`, `sin1`, `difficult`, `difficult2`, and
  `double_sine`.
- Optional retention of the StoSOO search tree and full evaluation history.
- Deterministic random seeding helper for reproducible examples/tests.

`plotStoSOO` and R S4/environment/graphics infrastructure are intentionally
omitted.

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example stosoo_guirland
fpm run --example poo_quadratic
```

A compiler-only validation script is also provided:

```text
scripts/test_gfortran.sh
```

or on Windows:

```text
scripts\test_gfortran.bat
```

## API sketch

```fortran
use oor

type(stosoo_options) :: options
type(stosoo_result)  :: result
real(dp) :: lower(2), upper(2)

lower = 0.0_dp
upper = 1.0_dp
options%stochastic = .false.
call stosoo(objective, lower, upper, 500, result, options)
```

The objective supplied to `stosoo` has the interface

```fortran
function objective(x) result(value)
   use oor, only : dp
   real(dp), intent(in) :: x(:)
   real(dp) :: value
end function objective
```

For `poo`, the objective receives a scalar `real(dp)` and is maximized on
`[0,1]`, matching the R package.

## License and provenance

The upstream package's `DESCRIPTION` declares `License: LGPL` without a
version qualifier.  This translation preserves that declaration.  The
complete supplied R package is retained under `original/OOR-master/`.

See `TRANSLATION_COVERAGE.md` and `UPSTREAM_PROVENANCE.md` for details.
