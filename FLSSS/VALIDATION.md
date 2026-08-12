# Validation

Validated with GNU Fortran 14.2.0.

## Strict serial build

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

All eight regression executables pass:

* `test_subset_search`
* `test_knapsack`
* `test_gap`
* `test_arbitrary`
* `test_controls`
* `test_v020_acceleration`
* `test_v030_fastpaths`
* `test_v040_mpat`

`test_v040_mpat` compares complete mPAT solution sets with the v0.2 int64 reference, exercises the public `engine='mpat'` path, checks automatic mPAT selection on a central comonotonic target, compares imposed-bound searches, and runs a set of deterministic randomized equivalence cases. It also checks that the triangular cache, bound-update and split diagnostics are active.

The v0.3 regression remains in the suite and continues to compare packed/PAT paths with the older reference engines and to verify decomposition behavior.

## Optimized warning-clean build

The complete library and all eight tests also pass with:

```text
-std=f2018 -O2
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

As part of v0.4, two remaining v0.3 expressions that combined an index-range test with an impure bound-test function using `.and.` were rewritten as nested tests. This avoids GNU Fortran's `-Wfunction-elimination` diagnostic and removes any dependence on short-circuit evaluation.

## OpenMP compatibility

The v0.3 decomposition/OpenMP layer is unchanged. The fast-path/decomposition regression also passes when the package is compiled with `-fopenmp`. mPAT itself is a serial search kernel in v0.4; decomposed parallel search continues to use the v0.3 packed exact engine.

## v0.4 mPAT performance

All measurements below use GNU Fortran 14.2.0 `-O2` in this environment. They are workload demonstrations, not universal guarantees.

### Central 24-dimensional search

`example/benchmark_mpat.f90`: 45 rows, 24 dimensions, subset size 10, exact central target.

```text
packed DFS: ~0.64 s, 4,168,901 states
v0.3 PAT:   ~0.059 s,   20,274 states
v0.4 mPAT:  ~0.045 s,   20,274 states

triangular entries: 405
triangular lookups:  708,462
packed/mPAT ratio:   ~14x
PAT/mPAT ratio:       ~1.3x
```

The PAT and mPAT state counts are identical. The v0.4 gain over PAT therefore comes from the triangular cached replacement sums and incremental `sumLB`/`sumUB` state, not from changing the represented search space.

### Larger central interval-state comparison

A separate 60-row, 16-dimensional, subset-size-10 central target produced approximately:

```text
v0.3 PAT:  14.31 s, 5,337,538 states
v0.4 mPAT: 10.55 s, 5,337,538 states
ratio:      1.36x
```

This confirms that the triangular cache continues to help when the interval-state tree is large.

## Earlier accelerators retained

The v0.3 broadword packed engine, v0.2 completion envelopes, arbitrary-precision k-sum lookup, and OpenMP decomposition remain present and regression-tested. The `dfs`, `packed`, and `pat` engines are retained for direct A/B comparisons with `mpat`.

## Fresh-archive check

The final release archive is extracted into a fresh directory, rebuilt from the archived sources under the strict serial flags above, and all eight regression executables are rerun before release.

## Source-format checks

No `src/`, `test/`, or `example/` Fortran source line exceeds 132 characters. No nonstandard free-form line-length option is required.
