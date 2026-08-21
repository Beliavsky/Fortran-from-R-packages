# Upstream dependency assessment

`trawl` 0.2.2 imports `DEoptim`, `rootSolve`, `Runuran`, `TSA`, and several
plotting packages. Only the first four participate in non-plotting computation.

## DEoptim

Version 0.2.0 embeds the supplied `DEoptim-fortran` 0.1.0 translation of
DEoptim 2.2-8. The translated DEoptim numerical engine is compiled directly
into the FPM source tree and is now used by:

- `fit_supigtrawl`
- `fit_lmtrawl`
- `fit_dexptrawl`

The wrapper uses the same controls selected by upstream `trawl`:
`strategy=2`, `NP=10*npar`, `CR=0.5`, `F=0.8`, `itermax=1000` by default,
and no tracing. The complete supplied dependency translation is retained under
`vendor/DEoptim-fortran-v0.1.0/` for provenance and licensing.

DEoptim in R draws from R's global RNG. The standalone Fortran DEoptim port has
its own RNG, so the trawl wrapper seeds that engine from one draw of the trawl
RNG. Consequently `set_trawl_seed()` makes DEoptim-backed trawl fits
reproducible, although streams are not seed-for-seed identical to R.

## rootSolve

A full `rootSolve` dependency is not required. `trawl` uses only
`rootSolve::uniroot.all` for scalar crossing detection. This is translated
locally using interval scanning plus bisection.

## Runuran

A full `Runuran` dependency is not required. `trawl` uses only the logarithmic
series RNG, which is implemented locally with the same distribution
parameterization.

## TSA

A full TSA dependency is not required for this package. `trawl` uses only an
univariate ACF call in its fitting routines; the corresponding demeaned sample
ACF calculation is implemented locally.

## Plotting imports

`MASS`, `squash`, `ggplot2`, `ggpubr`, and `graphics` are used by plotting code
or plotting branches and are intentionally omitted from the Fortran build.
