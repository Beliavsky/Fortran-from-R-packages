# Licensing

## gkwdist-derived Fortran source

The translated source under `src/` is derived from `gkwdist` 1.1.4 and is
provided under the MIT license, matching the upstream package declaration
`MIT + file LICENSE`. The upstream package is retained under `upstream/`.

## numDeriv-fortran

The supplied `numDeriv-fortran` translation is retained under
`vendor/numderiv-fortran/` under GPL-2.0-or-later. It is used by the derivative
regression test to independently verify the analytical derivatives. The
`gkwdist-fortran` library source does not import or call `numderiv`; therefore
its production library remains MIT-licensed. Programs that directly link and
use the GPL dependency must comply with that dependency's license.
