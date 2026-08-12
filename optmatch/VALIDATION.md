# Validation

Validation was performed with GNU Fortran 14.2.0 using strict runtime and
compile-time checks:

```text
-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
```

The environment used for this translation did not provide the `fpm` executable,
so the same FPM source/test tree was compiled directly with `gfortran` in module
dependency order.

## Deterministic suite

`test/test_optmatch.f90` exercises:

- optimal 1:1 matching;
- ordinary full matching;
- one-to-many and many-to-one full matching;
- exact restrictions and calipers;
- dense/sparse round trips;
- score, Euclidean, Mahalanobis, and rank-Mahalanobis distances;
- standardization scale;
- effective sample size;
- caliper sizing;
- `minControlsCap` / `maxControlsCap` style searches;
- integer utility translation.

The strict checked build passes the suite.

## Randomized solver validation

`test/test_random_fmatch.f90` generates 400 small random bipartite problems with
varying ratio constraints, forbidden edges, and selected omission fractions.
For each case it enumerates every possible subset of eligible match arcs and
computes the exact minimum objective satisfying the same network constraints.
The Fortran min-cost-flow result matched brute force on feasibility and objective
in all 400 cases.

This is particularly useful because it validates both the min-cost-flow kernel
and the nontrivial End/Sink capacity construction used by optimal full matching.
