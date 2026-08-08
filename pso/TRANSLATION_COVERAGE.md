# Translation coverage

## Translated computational code

The numerical behavior in `R/psoptim.R` is translated, including:

- SPSO2007 initialization and movement
- SPSO2011 initialization and rotationally structured movement
- informant-link generation and best-informant selection
- random or sequential asynchronous particle processing
- synchronous/vectorized particle processing
- inertia interpolation from `w0` to `w1`
- personal and global best tracking
- bound handling and zeroing velocity components that hit bounds
- Euclidean velocity clamping
- restart based on swarm diameter
- maximal restart and stagnation termination
- function-evaluation and iteration limits
- `fnscale` objective scaling/maximization
- `hybrid="on"` and `hybrid="improved"` behavior
- trace-statistics collection

The computational parts of `R/test.problem.R` and `R/test.result.R` are also
translated:

- all five benchmark objectives and gradients
- benchmark metadata/default bounds
- repeated PSO trials
- success-rate curves
- efficiency calculation
- objective summary statistics

## Intentionally omitted

The following R-specific or graphical code is not translated:

- S4 classes/method dispatch machinery
- `show()` formatting
- `plot()`
- `lines()`
- `points()`
- R namespace/registration infrastructure

The data represented by the S4 classes is represented by ordinary Fortran
derived types instead.

## Differences from the R runtime

### Hybrid local optimizer

Upstream `pso` delegates hybrid refinement to R's external
`optim(..., method="L-BFGS-B")`; that implementation is not part of the `pso`
package source.  To keep this FPM package standalone, the translation includes
a bounded limited-memory BFGS local optimizer with projected active bounds and
Armijo backtracking.  Consequently hybrid trajectories need not match R's
L-BFGS-B trajectories exactly.

When no analytical gradient is supplied, the Fortran local solver uses central
finite differences.  As documented by upstream for `maxf`, these numerical
finite-difference calls are not added to the PSO function-evaluation count.

### Random-number stream

The translation uses Fortran's intrinsic `random_number`.  It preserves the
same random constructions used by the R source, but it does not reproduce R's
RNG stream bit-for-bit.  Exact stochastic trajectories therefore differ even
when equivalent integer seeds are used.

### Timing

R's `proc.time()` returns user, system, and elapsed components.  The portable
Fortran benchmark summary stores CPU time and wall-clock elapsed time; it does
not attempt to synthesize a separate system-CPU component.

### Source behaviors preserved literally

Two potentially surprising upstream details are intentionally retained:

1. The helper named `rsphere.unif` does **not** draw a Gaussian/isotropic
   direction.  It draws positive `runif` components, normalizes that vector,
   and multiplies by a uniformly sampled radius.  The Fortran translation does
   the same.
2. In the upstream *vectorized SPSO2007* branch, the global-informant random
   coefficient uses `c.p`, not `c.g`.  The translation preserves that source
   behavior rather than silently correcting it.

The shifted Rosenbrock gradient is also translated literally from upstream,
including its `2*t0` term.  This is retained for source fidelity even though
it does not exactly differentiate the objective's `sum(x^2)` term.

The Ackley gradient adds a zero-norm guard so evaluating the gradient at the
exact optimum returns zero instead of forming the upstream expression's `0/0`.
