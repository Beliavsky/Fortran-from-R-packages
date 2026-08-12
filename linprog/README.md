# linprog-fortran 0.1.0

Modern Fortran/FPM translation of the computational code in R package
`linprog` 0.9-6 by Arne Henningsen.

## What is translated

- `solveLP()`'s original two-phase tableau simplex algorithm.
- Phase-I artificial-variable construction and Phase-II optimization.
- The package's pivot-selection heuristic.
- Basic-variable and full-variable results.
- Shadow prices / dual values.
- Objective-coefficient stability ranges and marginal-validity ranges.
- Constraint slack/free values.
- Optional explicit solution of the dual problem.
- The optional `lpSolve=TRUE` solver path, using the supplied modern Fortran
  translation of `lpSolve` as an FPM path dependency.
- `readMps()` and `writeMps()` including the same restricted MPS feature set as
  the R package (`L`/`G` rows and `UP` bounds; equality rows and `LO`/`FX`/`FR`
  bounds remain unsupported by the original parser).

R S3 printing and summary methods are presentation code and are not translated.

## Public API

```fortran
use linprog

type(linprog_control) :: control
type(linprog_result)  :: result

control%maximum = .true.
call solveLP(cvec, bvec, amat, result, control)
```

Constraint directions can be supplied as strings matching the R package:

```fortran
character(len=2) :: direction(3)
direction = ['>=', '>=', '<=']
call solveLP(cvec, bvec, amat, result, control, direction)
```

For integer direction constants, use `solveLP_dirs()` with `LINPROG_LE`,
`LINPROG_EQ`, and `LINPROG_GE`.

To use the supplied lpSolve translation:

```fortran
control%use_lpsolve = .true.
```

To solve the explicit dual problem as well:

```fortran
control%solve_dual = .true.
```

MPS I/O:

```fortran
type(mps_model) :: model
call writeMps('problem.mps', cvec, bvec, amat, 'Example')
call readMps('problem.mps', model, solve=.true., maximum=.false.)
```

## Important compatibility behavior

The original internal `solveLP()` implementation is documented by its author as
unsafe for equality constraints.  It maps `=`/`==` to a zero sign in its tableau,
which effectively removes the equality coefficients.  This behavior is
preserved in the internal Fortran path for source compatibility.

For equality constraints, use:

```fortran
control%use_lpsolve = .true.
```

The regression suite reproduces the original package example where the internal
path returns objective 1998 while the lpSolve-backed path correctly returns 1404.

The original package also contains two literal sensitivity placeholders (`99`
and `77`) for nonbasic decision variables in minimization problems.  They are
preserved rather than silently replaced with invented sensitivity values.

## FPM dependency

The supplied dependency is vendored under `vendor/lpsolve-fortran` and declared
in `fpm.toml` as a local path dependency.  No external R installation is needed.

## Validation

Five regression executables cover:

1. Both published `linprog` examples, including Phase I and detailed dual data.
2. Internal-vs-lpSolve equality behavior.
3. Explicit dual solving with both solver backends.
4. `<=`, `>=`, and equality direction handling.
5. MPS write/read round-trip plus `G` rows and `UP` bounds.

See `VALIDATION.md` for build details.
