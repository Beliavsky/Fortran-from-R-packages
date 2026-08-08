# Translation coverage

## Upstream computational code translated

The upstream package contains one computational R source file, `R/soma.R`.
The following functionality is translated:

- `bounds`
- `all2one`
- `t3a`
- `pareto`
- population initialization from bounds
- user-supplied initial populations
- initial cost evaluation
- All To One leader selection and migration
- T3A leader/migrant pool sampling
- T3A adaptive perturbation probability
- T3A adaptive step lengths
- Pareto cost ranking
- Pareto leader and migrant percentile pools
- Pareto sinusoidal perturbation probability
- Pareto sinusoidal step lengths
- per-dimension perturbation masks
- candidate points along migration trajectories
- overshoot beyond the leader
- random replacement of out-of-bounds coordinates
- selection of each migrant's lowest-cost candidate
- absolute and relative population-cost separation termination
- migration limits
- no-op migration behavior
- leader-cost history
- evaluation-count history

The Fortran code preserves the R implementation's column-major population
layout: parameters are rows and individuals are columns.

## Source behaviors deliberately preserved

The R implementation generates a full
`n_parameters * population_size * n_steps` random repair array for every
counted migration even when only a subset of individuals migrates. The
Fortran translation generates the same number of repair random values and
uses the corresponding linear-index values for out-of-bounds candidates.

Pareto's perturbation mask is also preserved exactly in form: dimensions not
selected by the Bernoulli perturbation use `progress` rather than zero, so
later Pareto migrations still move those dimensions more slowly.

The evaluation-history convention is preserved: the history starts at zero
cost-function evaluations even though the initial population has already been
evaluated, because this is how the R object used by `plot.soma` is constructed.

## Omitted non-computational code

- `plot.soma`
- plotting code in the demo and README
- `reportr` logging/output
- S3 class metadata and R object plumbing
- roxygen/man-page generation

## Differences from the R runtime

### Random-number stream

Fortran uses `random_number` and a Fortran sampling-without-replacement
routine. Therefore `set.seed(18)` in R and `call soma_set_seed(18)` in Fortran
do not produce identical trajectories. The algorithmic distributions and
random draw structure are retained where practical.

### Bounds diagnostics

The R `bounds()` helper reports some malformed bounds through `reportr`
warning-level assertions. The Fortran optimizer reports a nonzero `status`
for a lower bound greater than its upper bound because there is no R warning
condition system.

### Extra objective arguments

R's `soma(..., ...)` can forward arbitrary extra arguments to the objective.
Fortran callbacks instead use normal Fortran mechanisms such as module state
or derived types.

### Reporting

Progress and informational console messages are omitted. They do not affect
the optimizer state.
