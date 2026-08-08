# Translation coverage

Upstream: `calibrar` 0.9.0.

## Directly translated calibrar algorithms

### `R/calibrar-gradient.R`

Translated forward, backward, central and Richardson finite-difference gradients, including calibrar's default step-size rules and the near-zero derivative shortcut in Richardson extrapolation. R parallel execution is omitted; serial and parallel R paths calculate the same mathematical estimates.

### `R/calibrar-AHR_ES.R`

Translated the AHR-ES numerical engine:

- truncated/clipped Gaussian population creation;
- mean individual inserted as population member 1;
- global and per-objective ranking;
- log recombination weights;
- hierarchical per-objective mean/dispersion smoothing;
- CV-based objective recombination weights;
- evolution paths `pc` and `ps`;
- diagonal variance adaptation;
- global step-size adaptation;
- multi-objective weighted aggregation;
- termination rules 0--4 and trace arrays.

Restart RDS files and progress-printing/timing infrastructure are omitted.

### `R/calibrar-core-calibrar.R` and `R/calibrar-core-optim2.R`

Translated the numerical orchestration concepts:

- active parameter masks;
- sequential phases;
- bound checking/clamping;
- per-phase stochastic replicate averaging;
- scalar and multi-objective calibration;
- common optimizer result and options types.

R's `relist`/list-shaped parameters are represented as flat real arrays. This is the natural representation for numerical Fortran APIs.

### `R/calibrar-fitness*.R` and objective construction

Translated the package-owned fitness/error formulas and weighted aggregation. `calibration_objFn`'s numerical data comparison is represented by `calibration_term` plus flat observed/simulated arrays rather than an R closure over named lists/data frames.

### `R/calibrar-random.R`

Translated normal, truncated-normal, multivariate-normal density and Gaussian-kernel functionality. `rtnorm_matrix` preserves calibrar's `rtnorm2` policy: dimensions with small excluded Gaussian mass use a normal draw clipped to bounds; more strongly truncated dimensions use rejection/tail sampling.

### `R/calibrar-splines.R`

Translated spline parameter interpolation and Gaussian-kernel computation. Plotting is omitted.

### `R/calibrar-stopping.R`

Translated `smooth_stop2`, `smooth_stop3`, `smooth_stop4`, and `N_stop` logic.

## External optimizer wrappers

Much of `calibrar-wrapper*.R` only adapts arguments/results to algorithms implemented in other R packages. Those external implementations are not part of the calibrar source archive.

The standalone Fortran dispatcher therefore distinguishes package-owned algorithms from compatibility implementations:

| R method | R implementation lives in | Fortran standalone behavior |
|---|---|---|
| AHR-ES | calibrar | Direct translation |
| BFGS | stats | Self-contained BFGS analogue |
| CG / Rcgmin | stats / optimx | Self-contained nonlinear CG analogue |
| L-BFGS-B / LBFGSB3 / nlminb | stats / lbfgsb3c / stats | Projected BFGS compatibility path |
| Nelder-Mead / nmk / nmkb | stats / dfoptim | Self-contained Nelder-Mead path |
| hjn / hjk / hjkb / mads | optimx / dfoptim | Hooke-Jeeves compatibility path |
| spg | BB | Self-contained spectral projected gradient path |
| SANN / genSA | stats / GenSA | Self-contained simulated-annealing path |
| CMA-ES, DE, soma, genoud, PSO, hybridPSO | external packages | Integration point; not falsely substituted |

Exact behavior for an external optimizer requires linking an appropriate Fortran translation of that external package. This is intentionally different from copying or rebranding another method as the requested algorithm.

## R-specific / non-computational code omitted

- S3 print/summary/plot/predict methods and plotting code;
- `foreach`, PSOCK/FORK cluster setup, and parallel filesystem work directories;
- RDS restart/result serialization;
- command-line/configuration convenience parsing;
- R data-frame/list class machinery;
- demo-file generation and package-vignette infrastructure.

The complete originals remain under `original/calibrar-master/`.

## Known numerical differences

- External R optimizer implementations are not bit-for-bit reproduced by the standalone compatibility methods listed above.
- Fortran RNG streams are not R RNG streams, so stochastic trajectories differ even when the same integer seed is supplied.
- `gaussian_kernel_2d` uses an evenly spaced requested grid rather than R's `pretty()` expansion.
- The spline implementation is a standalone cubic not-a-knot solver rather than calling R's `splinefun` implementation.
- R's disk/parallel replicate evaluation order is not reproduced; replicate averaging is preserved numerically in serial execution.
