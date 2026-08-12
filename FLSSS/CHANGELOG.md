# Changelog

## 0.4.0

* Added `search_md_i8_mpat`, a native mPAT-style engine for comonotonic integerized multidimensional search. It keeps lower/upper index vectors and cached packed `sumLB`/`sumUB`, alternates bound tightening to a fixed point, and branches on the smallest nonzero index gap.
* Added a packed triangular consecutive-sum cache analogous to the original FLSSS `TriM/M` data. Every consecutive row sum for lengths `1..len` is precomputed once, allowing lower/upper candidate tests to replace a whole consecutive prefix/suffix in packed-lane time rather than resumming all selected rows.
* Stored triangular entries as `value(:,entry)` so all lanes used by a lookup are contiguous in Fortran memory.
* Child mPAT states carry `sumLB` and `sumUB` incrementally across splits; bound changes update these sums by subtracting/adding only changed rows.
* Added public `engine="mpat"`. `engine="pat"`, `"packed"`, and `"dfs"` remain available as v0.3/v0.2 reference backends.
* Broadened `engine="auto"` so moderate/large comonotonic central-target problems can select mPAT instead of limiting interval-state search to near-boundary targets.
* Added diagnostics `tri_entries`, `tri_lookups`, `bound_updates`, and `mpat_splits`.
* Added `test_v040_mpat` with complete-set, randomized, public-API, central-target auto-selection, and imposed-bound comparisons against the v0.2 reference.
* Added `benchmark_mpat`, comparing packed DFS, the v0.3 PAT engine, and the v0.4 mPAT triangular-cache engine on the same central target.
* Removed two remaining `.and.` calls around impure bound-test functions in the v0.3 PAT code; explicit nested tests avoid compiler-dependent short-circuit/function-elimination diagnostics at `-O2 -Werror`.

## 0.3.0

* Added an exact broadword packing layer for native `integer(int64)` multidimensional search. Multiple dimensions share 64-bit lanes with per-field guard bits; lower/upper interval comparisons are evaluated fieldwise without cross-field carries. Full exact arithmetic is retained.
* Added `search_md_i8_packed`. The high-level integerized API now accepts `engine="auto"`, `"packed"`, `"pat"`, or `"dfs"`; `dfs` is the unchanged v0.2 reference engine.
* Added a PAT-inspired comonotonic interval-state engine. It maintains lower/upper index vectors, repeatedly tightens them with monotone lower/upper feasibility searches, branches on the position with the smallest remaining gap, and uses the packed comparison layer for state tests.
* Added a boundary-aware `auto` heuristic: packed DFS is the normal choice; PAT is selected only for comonotonic integerized problems whose target lies close to a feasible edge, where interval tightening is most effective.
* Added safe fixed-cardinality zero-minimum shifting before packed search. Targets are adjusted exactly in int64 arithmetic; unsafe-overflow cases fall back to the v0.2 engine.
* Added optional OpenMP execution for multidimensional and arbitrary-precision decomposition objects. Without `-fopenmp`, the same source is a serial build.
* `maxCore` compatibility wrappers now use decomposition/parallel execution when and only when the package was compiled with OpenMP. Ordinary builds retain the serial behavior rather than paying decomposition overhead.
* Added `mflsss_par_integerized_parallel`, `mflsss_decomp_run`, and `arb_flsss_decomp_run`.
* Added search diagnostics `bound_states`, `partitions_run`, `packed_lanes`, and `engine`.
* Serial decomposition runners now enforce a shared wall-clock budget across partitions.
* Added `test_v030_fastpaths` covering packed-vs-v0.2 solution equivalence, PAT equivalence, public engine selection, automatic PAT selection near a boundary, and serial/parallel decomposition equivalence.
* Added `benchmark_packed`, `benchmark_pat`, and `benchmark_parallel` examples.

## 0.2.0

* Added exact completion-envelope pruning for 1-D and multidimensional subset-sum search. The envelopes precompute the minimum and maximum feasible completion from every `(position,start)` state and work for unsorted data and position-specific index bounds.
* Added a native `integer(int64)` multidimensional mining path. Integerized searches no longer convert the integerized matrix, target and tolerance back through `real(dp)` before search.
* Turned `ksumHash` into an actual arbitrary-precision lookup accelerator: generated tables contain sorted 64-bit tuple fingerprints, while lookup verifies the full canonical decimal-integer tuple to make collisions harmless.
* `arb_flsss(..., given_ksum=...)` now uses the supplied k-sum table by splitting every candidate uniquely into a left prefix and the final `k` indices.
* Fixed a v0.1.0 arbitrary-precision DFS state bug: a host-shared temporary sum buffer could contaminate sibling recursive branches. Recursive partial sums are now saved in invocation-local storage.
* `decompose_mflsss` and `decompose_arb_flsss` now honor `approx_ninstance` by creating exact contiguous first-index ranges. Resume objects cover the complete search space without overlap.
* Added search diagnostics `pruned`, `hash_lookups`, and `hash_candidates` to `subset_solutions`.
* Added `test_v020_acceleration` plus `benchmark_bounds` and `benchmark_ksum` examples.

## 0.1.0

* Initial modern Fortran/FPM translation of the exported FLSSS 9.2.8 computational surface.
* Added exact 1-D, multiset and multidimensional subset-sum search.
* Added bounded and variable-cardinality search plus conjugate 1-D search.
* Added decomposition/resume types.
* Ported active FLSSS integerization rules.
* Added multidimensional knapsack, 0/1 branch-and-bound and DP knapsack.
* Added exact maximum/minimum GAP and GA heuristic.
* Added arbitrary-precision decimal-string search and exact k-sum tables.
* Added R-name compatibility wrappers.
* Retained the complete original package under `original/`.
