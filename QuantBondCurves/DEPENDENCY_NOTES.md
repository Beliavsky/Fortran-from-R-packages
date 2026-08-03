# Dependency notes

The R package imports `lubridate`, `quantdates`, and `Rsolnp`.

## Date dependencies

The required date arithmetic and day-count calculations were translated into native Fortran. Location-specific holiday databases were not copied.

## Rsolnp attachment

The supplied `Rsolnp-fortran` archive declares `GPL-2.0-only`. The upstream `QuantBondCurves` package declares `GPL (>= 3)`, represented here as `GPL-3.0-or-later`.

GPL-2.0-only and GPL-3.0-or-later code cannot safely be distributed as one statically linked combined FPM executable. Consequently, this package does not copy or link the attached dependency. Curve and basis calibration use a clean native bounded Nelder-Mead optimizer under GPL-3.0-or-later.

This is a packaging-license decision, not a claim that the mathematical SOLNP algorithm cannot be implemented independently under another license.
