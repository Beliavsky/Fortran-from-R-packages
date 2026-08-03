# Algorithm

For visiting parameter `qv`, iteration age `t`, and initial temperature `T0`,
the translated temperature schedule is

```text
T(t) = T0 * (2^(qv-1) - 1) / ((t+1)^(qv-1) - 1).
```

When the temperature falls below `temp_restart`, the schedule age is restarted.
This retains GenSA's re-annealing behavior while keeping an explicit total
outer-iteration counter.

The visiting step follows the upstream `visita` implementation. It forms a
scale from `qv`, `T`, gamma functions, and two independent Gaussian deviates.
The resulting heavy-tailed proposal is wrapped periodically into each finite
box interval. During the first `dimension` moves of each default Markov chain,
all coordinates change. During the next `dimension` moves, one coordinate
changes at a time.

A downhill proposal is always accepted. For an uphill increase `delta`, the
acceptance probability is

```text
p = [1 + (qa - 1) * delta / Tqa]^(1 / (1 - qa)),
```

when the bracket is positive, and zero otherwise. `Tqa` is the current
visiting temperature divided by the current schedule age.

The best point is optionally polished by a bounded local search:

- smooth objectives use projected BFGS with bound-aware central differences
  and an Armijo line search;
- nonsmooth objectives use a bounded coordinate-pattern search with shrinking
  step lengths.

Nonfinite objective values are converted to a large finite sentinel and cannot
replace the incumbent. Every actual objective evaluation is counted.
