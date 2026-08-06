# Dependency and license notes

## Compiled dependencies

The source tree vendors the computational modules from these completed Fortran
translations:

- `ghyp-fortran`: generalized-hyperbolic mathematics used by
  `tsdistributions`
- `tsdistributions-fortran` 0.1.0: standardized innovation distributions,
  fitting support, numerical optimization, and moments
- `tsgarch-fortran` 0.1.0: univariate GARCH models, filtering, estimation,
  simulation, and forecasting

These sources are compatible with the GPL-2.0-only license selected for this
combined package.  Their original archives are retained under
`provenance/dependencies`.

## User-supplied nloptr translation

`provenance/dependencies/nloptr-fortran(2).zip` is preserved unchanged for
reference.  It is not included in `scripts/common.sh`, the FPM source graph, or
any `use` statement.

The supplied translation declares LGPL-3.0-or-later.  GPL version 2 only is not
compatible with LGPL version 3 for static combination, so linking it into this
GPL-2-only port would create a license conflict.  The bounded optimizer already
available in the GPL-2-compatible `tsdistributions` translation is used instead.

## Replaced R dependencies

- `Rsolnp`, `nloptr`: bounded derivative-free optimizer from
  `tsdistributions-fortran`
- `numDeriv`, `sandwich`: finite-difference derivatives and covariance
  calculations
- `RcppArmadillo`, `RcppParallel`: native Fortran arrays and serial numerical
  kernels
- `RcppBessel`: generalized-hyperbolic special functions from `ghyp-fortran`
- `future`, `future.apply`: deterministic serial loops
- `xts`, `zoo`, `data.table`: plain numeric arrays; date indexes are left to
  callers

No code from the provenance-only `nloptr` archive is compiled.
