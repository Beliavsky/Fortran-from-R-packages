# Notices and attribution

This is a computational modern-Fortran translation of SpatialExtremes 2.1-0.

Upstream package:

- Author: Mathieu Ribatet
- Contributor: Richard Singleton
- Contributor: R Core team
- License: GPL (>= 2)
- Upstream package copyright/license text is retained verbatim in `LICENSE-SpatialExtremes.txt` and `upstream/SpatialExtremes-Copyright`.

The upstream package cites and implements methods for spatial extremes including work by Schlather, Smith, Brown and Resnick, Padoan/Ribatet/Sisson, Davison/Padoan/Ribatet, Dombry/Engelke/Oesting and related authors. The original R/C sources are retained under `upstream/` to preserve detailed algorithm provenance. Version 0.2.0 adds Fortran translations derived in particular from `standardErrors.c`, `maxLinear.c`, `randomlines.c`, `turningbands.c`, the TBM sections of `simschlather.c`/`simgeometric.c`/`simextremalt.c`, `condsimMaxStab.c` partition logic, and `latentVariable.c` DIC/latent likelihood kernels. Version 0.3.0/0.4.0 additionally derives numerical translations from `latentVariable.c`/`mcmc.c`, `circulant.c`, `condsimMaxStab.c`, `newCode.c`, `utils.c` (`gev2frechTrend`), and the model-specific standard-error paths in `standardErrors.c`/`standardErrorsCommonPart.c`.

`src/r_mod.F90` is based on the separately supplied MIT-licensed `r_mod.f90`; its original is retained as `upstream/r_mod-original.f90` and remains under those MIT terms.

The upstream file `src/fft.c` is copied from R and carries its own embedded R copyright/license statements.  It is retained only in the archival upstream source tree; the Fortran public API in this release does not compile or redistribute that C implementation as translated code.
