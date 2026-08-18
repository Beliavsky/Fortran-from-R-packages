# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `mcmc` 0.9-8.
- Ported random-walk Metropolis with all three scale modes, batch acceptance,
  output callbacks, and debug trajectories.
- Ported serial and parallel tempering, including neighbor Hastings correction
  and swap acceptance accounting.
- Ported morphometric radial transformations and morphed Metropolis.
- Ported `initseq` PAVA/convex-sequence variance estimators.
- Ported `olbm` overlapping batch means covariance.
- Added strict behavioral tests and standalone examples.
