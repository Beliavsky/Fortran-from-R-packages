# calibrar-fortran

Modern Fortran translation of the computational core of the R package `calibrar` 0.9.0 (Ricardo Oliveros-Ramos), packaged for the Fortran Package Manager (FPM).

The translation focuses on numerical and calibration functionality. R S3 classes, plotting, `foreach`/cluster orchestration, RDS serialization, command-line helpers, and data-frame/file-management plumbing are intentionally not translated.

## Main functionality

- `calibrate` -- sequential/phased scalar calibration with active-parameter masks and per-phase stochastic replicates.
- `calibrate_multi` -- phased multi-objective calibration through AHR-ES.
- `ahres` / `ahres_scalar` -- Adaptive Hierarchical Recombination Evolutionary Strategy.
- `optim2` / `optimh` -- common Fortran optimization interface.
- Numerical gradients: forward, backward, central, Richardson.
- Calibration fitness functions: `norm2`, `lnorm2`, `lnorm3`, `lnorm4`, `lnorm4b`, Poisson, penalties, and multinomial-style error.
- Flat-array calibration objective construction through `calibration_term`.
- Truncated-normal generation following calibrar's clipped/exact hybrid rule.
- Multivariate-normal density and 2-D Gaussian-kernel grid.
- Cubic spline parameter interpolation (`spline_par`).
- calibrar stopping criteria (`smooth_stop2`, `smooth_stop3`, `smooth_stop4`, `n_stop`).

## Build

```sh
fpm build
fpm test
```

Examples:

```sh
fpm run --example basic_optim
fpm run --example ahres_multiobjective
```

## Optimizer wrappers

The R package delegates many methods to other R packages. Their implementations are not source code belonging to calibrar. This translation does not claim those external algorithms are exact calibrar translations.

The standalone `optim2` dispatcher provides native Fortran implementations/analogues for BFGS/projected BFGS, nonlinear CG, Nelder-Mead, Hooke-Jeeves, spectral projected gradient, simulated annealing, and calibrar's own AHR-ES. Compatibility aliases are documented in `TRANSLATION_COVERAGE.md`.

External package methods such as CMA-ES, DEoptim, pso, soma, rgenoud, and exact minqa/dfoptim/optimx variants should be connected through separate Fortran translations/adapters when exact external-package behavior is required.

## Licensing

The upstream package declares `License: GPL-2`. This translation is distributed under GPL-2.0-only, and the complete supplied upstream package is retained in `original/calibrar-master/` for provenance.
