# Algorithm notes

## Pareto and hypervolume

GPareto's C++ dominance routines are translated directly in semantics. The
Fortran hypervolume implementation is an exact recursive sweep for minimization
problems and removes the runtime dependency on `emoa`.

## Expected hypervolume improvement

For two objectives, `ehi_2d_values` follows the upstream
`EHI_2d_wrap_Rcpp` cell decomposition, including the same `exipsi` marginal
integrals and dominated-cell correction. Zero predictive variance is handled
explicitly by returning deterministic hypervolume improvement.

For more objectives, `crit_ehi` uses sample-average approximation, as upstream
GPareto does. `crit_qehi` draws the whole candidate batch jointly within each
objective using DiceKriging's predictive covariance before computing the
hypervolume increment.

## Gaussian-process layer

The package vendors the GPL-3 option of the translated DiceKriging numerical
core. `gp_model` adds a trend-basis tag so GPareto criteria can construct the
correct new-data design matrix automatically. `gparetoptim` can call genuine
DiceKriging hyperparameter re-estimation after each new observation.

## SUR

Rather than translating the R list/precomputation structure and the specialized
`EEV.2D.computation` expressions literally, `crit_sur` performs the equivalent
one-step Gaussian conditioning directly. For each candidate and objective it
forms the joint posterior covariance between integration points and the
candidate, samples the future standardized observation, updates integration
means/variances analytically, recomputes probability of non-domination, and
averages the integrated Bernoulli uncertainty. The returned value is current
uncertainty minus expected future uncertainty.

## Optimization and designs

The R package delegates continuous criterion optimization to `rgenoud` or
`pso`; those external engines are not part of GPareto. This port supplies a
bounded differential-evolution optimizer and a discrete-search routine.

The integration-design API provides Halton low-discrepancy, Monte Carlo and SUR
importance modes. The name `sobol` is accepted as an alias for Halton so callers
can retain upstream control strings without pulling in randtoolbox's 1111-D
Sobol table.
