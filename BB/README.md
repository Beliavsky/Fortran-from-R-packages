# BB-fortran

Modern Fortran translation of the computational code in the R package
`BB` 2026.1.0, with an FPM package layout.

The upstream BB package is by Ravi Varadhan and Paul Gilbert.  This port keeps
its GPL-3 license and includes the original package snapshot under
`original/BB-master/` for provenance.

The R package's `projectLinear` routine depends on `quadprog`.  The three
library modules required from the user-supplied `quadprog-fortran` translation
are vendored in `src/` unchanged, and its GPL-2.0-or-later license is preserved
under `LICENSES/`.

## Implemented computational routines

- `spg`: spectral projected gradient optimization, methods 1, 2, and 3
- box-constrained SPG
- user-defined projection SPG
- linear equality/inequality constrained SPG through `quadprog`
- `sane`: spectral residual nonlinear-equation solver
- `dfsane`: derivative-free SANE nonlinear-equation solver
- `bboptim`: retry driver over BB step methods and nonmonotone-memory sizes
- `bbsolve`: retry driver over DF-SANE configurations
- `multistart_solve`
- `multistart_optimize`
- box, custom-projection, and linear-constraint multistart optimization
- optional Nelder-Mead startup for SANE/DF-SANE
- optional BFGS residual-minimization fallback for SANE/DF-SANE
- `project_box` and `project_linear`

Plotting, R object manipulation, R error-catching machinery, vignette
presentation code, and R-specific namespace/registration code are not ported.

## Build

```text
fpm build
fpm test
```

Run the examples with:

```text
fpm run --example spg_rosenbrock
fpm run --example solve_system
```

## Callback API

The public umbrella module is `bb`.

A scalar objective is a Fortran function:

```fortran
function objective(x) result(f)
  use bb, only: dp
  real(dp), intent(in) :: x(:)
  real(dp) :: f
end function objective
```

An analytic gradient is a subroutine:

```fortran
subroutine gradient(x, g)
  use bb, only: dp
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: g(:)
end subroutine gradient
```

A nonlinear system callback is:

```fortran
subroutine equations(x, f)
  use bb, only: dp
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: f(:)
end subroutine equations
```

A custom projection callback is:

```fortran
subroutine projection(x, projected, ok)
  use bb, only: dp
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: projected(:)
  logical, intent(out) :: ok
end subroutine projection
```

Host association can be used for the equivalent of R's `...` and
`projectArgs`.

## Main examples

Unconstrained optimization:

```fortran
use bb, only: dp, spg, spg_result

type(spg_result) :: fit
fit = spg(x0, objective, gr=gradient)
```

Box constraints:

```fortran
fit = spg_box(x0, objective, lower, upper, gr=gradient)
```

Linear constraints use rows of `A` and satisfy

```text
A * x >= b
```

with the first `meq` rows treated as equalities:

```fortran
fit = spg_linear(x0, objective, a, b, meq, gr=gradient)
```

This row-oriented `A` convention deliberately matches BB's `projectLinear`
API rather than the column-oriented internal convention of `quadprog`.

## Multistart orientation

R stores starting points by matrix rows.  The Fortran port stores each starting
point in a column, so `starts(:,k)` is starting point `k`.  The returned
`multistart_result%par(:,k)` has the same convention.

## Numerical notes

The core SPG, SANE, and DF-SANE iteration formulas and default tolerances follow
the R source.  The R package delegates optional Nelder-Mead and L-BFGS-B work
to `stats::optim`; this standalone port supplies local Nelder-Mead and dense
BFGS equivalents so it has no R dependency.  These fallback paths are therefore
algorithmically equivalent in purpose but not expected to reproduce R's exact
iteration counts.

See `TRANSLATION_COVERAGE.md` and `VALIDATION.md` for details.
