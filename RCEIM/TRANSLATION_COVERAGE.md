# Translation coverage and differences

## Directly translated computational behavior

`ceimOpt` is translated as `ceim_optimize` with the same main algorithm:

1. uniformly initialize a bounded population;
2. evaluate the objective and multiply by `mfactor` (`+1` minimize, `-1` maximize);
3. sort by ascending internal score;
4. retain `N_elite` points;
5. update elite best/mean/sample-sd histories;
6. update `mu` with `alpha` and `sigma` with the dynamic `alpha_d` expression;
7. sample the next population from independent normal marginals;
8. carry `N_super` elite members into the next generation;
9. optionally replace non-super members with the source's broad Gaussian "chaos" population;
10. clip all coordinates to their bounds;
11. stop on the source's iteration, fractional-sigma, no-improvement, or identical-elite-score conditions.

R's `stats::sd` behavior is matched with the sample standard deviation (`n-1` denominator).

## Deliberate non-R-runtime differences

### Objective callback

R resolves a function object/name dynamically with `match.fun`. Fortran uses an explicit procedure interface. This removes R's dynamic name lookup without changing the numerical algorithm.

### Random-number stream

The translation uses the compiler's `random_number` generator plus a Box-Muller normal transform. Fixed Fortran seeds are reproducible within a given compiler/runtime, but are not intended to reproduce R's `runif`/`rnorm` stream bit-for-bit.

### Parallel evaluation

The upstream optional `parallel::mclapply` path is R-runtime orchestration, not a separate numerical algorithm. The FPM library evaluates objectives serially so it has no OpenMP or process-runtime dependency. Callers can parallelize independent runs externally.

### Plotting and interaction

`plotConvergence`, `plotResultDistribution`, `plotEliteDistrib`, `overPlotErrorPolygon`, graphics-device calls, and `handIterative/readline` are omitted.

### Generic `sortDataFrame`

Only the behavior used computationally by `ceimOpt` is translated: stable ascending row sorting by the score column. A generic R data-frame sorting abstraction has no natural Fortran analogue and is unnecessary for the optimizer.

## Upstream source quirks handled defensively

The R source initializes `chaosCounter` but later uses `chaosGenCounter`. In ordinary runs `chaosGenCounter` is created when the first best score differs from zero; unusual first-generation cases can otherwise reference an undefined variable. The Fortran translation initializes the intended chaos counter explicitly.

`Nnew` is assigned in R only inside the sigma-update block but is later also used by chaos generation. The Fortran translation defines `Nnew = Ntot-N_super` at initialization.

R can produce a vector `maxSigIdx` when multiple sigmas tie, which can make `fracsig` non-scalar. The Fortran translation uses the first maximum, consistent with the intended scalar stopping statistic.

If the selected `mu` is exactly zero, R's `sig/abs(mu)` produces Inf (or NaN for 0/0). The translation uses a very large value for positive-sigma/zero-mu and zero for zero/zero, preventing a floating exception while preserving the intended "not yet fractionally converged" interpretation.

The R source technically permits `N_elite=1`, after which `sd()` is `NA` and the algorithm ceases to be numerically well-defined. The Fortran API requires at least two elite members.

The R loop performs no iteration when `epsilon=0` because `fracsig` is initialized to `10*epsilon`; it then reaches an undefined `elite` object on return. The Fortran API therefore requires `epsilon > 0` rather than returning an uninitialized result.

## Maximization result sign

Upstream returns the elite `S` value, where `S=-f(x)` for maximization. The Fortran result preserves that quantity as `result%score` but additionally reports the natural objective as `result%value`, which is usually what a Fortran caller expects.
