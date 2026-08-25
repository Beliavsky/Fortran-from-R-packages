# GPArotation-fortran

**Official CRAN title:** Gradient Projection Factor Rotation

Modern Fortran/FPM translation of the numerical core of the R package
`GPArotation` 2026.8-2.

The library implements gradient-projection factor rotation for orthogonal and
oblique transformations, the package's exported rotation criteria, random
starts, Lp rotation, EIV/echelon identification, normalization weights, and
non-graphical fit/simple-structure diagnostics.

## Main modules

- `gpa_api`: convenience API corresponding to the exported R rotation wrappers.
- `gpa_rotation`: `gpforth`, `gpfoblq`, random-start and Lp engines.
- `gpa_criteria`: objective functions and analytical gradients (`vgQ.*`).
- `gpa_transforms`: EIV and echelon rotations.
- `gpa_diagnostics`: residuals, SRMR/RMSEA and simplicity indices.
- `gpa_linalg`: self-contained matrix utilities used by the algorithms.

All floating-point calculations use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

and all project source is free-form Fortran 2018.

## Example

```fortran
use gpa_kinds, only: dp
use gpa_rotation
use gpa_criteria, only: criterion_options

type(rotation_result) :: fit
type(rotation_options) :: opt
type(criterion_options) :: crit
real(dp) :: a(6,2)

! ... set a ...
opt%eps = 1.0e-7_dp
call gpforth(a, "varimax", fit, crit, opt)
```

The higher-level `gpa_api` module also provides procedures such as `varimax`,
`quartimin`, `geomin_t`, `geomin_q`, `bifactor_t`, `bifactor_q`, `equamax`,
`parsimax`, `target_t`, `target_q`, `lp_t`, and `lp_q`.

## Build

With FPM:

```text
fpm test
fpm run --example basic_rotation
```

The translation has also been validated directly with GNU Fortran 14 using:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all
```

## Scope

Plotting (`plot.GPArotation`, trajectory/landscape plots and
`plot2fOrthComparison`) and R-specific S3/factanal/model-object presentation
are not reproduced. Their underlying matrix calculations are represented where
they are useful independently of R.

See `API_MAP.md` and `PORTING_NOTES.md` for details.
