# Porting notes

## Scope

All exported numerical routines in the original package were translated.
R package infrastructure, roxygen-generated help machinery, pkgdown setup, and
R testthat plumbing are retained only in `original/` and are not compiled.

## Numerical replacements

The R implementation uses `statmod::gauss.quad` in its implied-volatility
workflow. The Fortran port instead uses self-contained bracket expansion and
bisection. CEV evaluation uses internal regularized-gamma and noncentral
chi-square routines.

## Source mapping

- `R/bachelier_price.R`, `R/bachelier_impvol.R` -> `src/fer_vanilla.f90`
- `R/blackscholes_price.R`, `R/blackscholes_impvol.R` -> `src/fer_vanilla.f90`
- `R/cev.R` -> `src/fer_cev.f90`, `src/fer_special.f90`
- `R/sabr.R` -> `src/fer_sabr.f90`
- `R/spread.R` -> `src/fer_spread.f90`
- Public aggregation -> `src/fer.f90`

## Licensing

Although the upstream repository includes a GPL-3 text, its package DESCRIPTION
states `GPL (>= 2)`. The port therefore uses the SPDX identifier
`GPL-2.0-or-later` and includes complete GPL-2.0 and GPL-3.0 texts.
