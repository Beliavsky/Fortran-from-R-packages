# Translation notes

## Source

* R package: FLSSS 9.2.8
* Package license: GPL-3
* Native dependencies in R package: Rcpp and RcppParallel
* Bundled third-party source: xxHash (`xxhash.h`, BSD-2-Clause)

## Search architecture through v0.4

The original package contains several generations of highly templated C++ exact-search engines. The Fortran project translates the exported mathematical problems into explicit numerical kernels and keeps successive engines available for comparison.

### v0.2 reference engines

`search_1d`, `search_md`, and `search_md_i8` use exact dynamic completion envelopes. For every remaining selection position and possible starting index, the code stores minimum and maximum achievable continuation sums. A state is rejected as soon as any constrained target interval cannot intersect its achievable range. These routines work for arbitrary row order and position-specific lower/upper bounds.

`search_md_i8` remains available through `engine='dfs'` and is the numerical reference for the v0.3/v0.4 integerized engines.

### v0.3 packed int64 engine

`flsss_packed` builds an exact broadword plan for nonnegative fixed-cardinality data.

For dimension `j`, the plan reserves enough data bits to hold `len*max(v(:,j))` and one additional guard bit. Fields are greedily grouped into positive 64-bit lanes while leaving the sign region unused. Because a complete subset sum is guaranteed to fit below the guard bit, adding packed row values cannot carry into an adjacent field.

Fieldwise comparisons use the standard guard-bit subtraction idea:

* `(sum | guards) - lower` exposes a guard bit for every field satisfying `sum >= lower`;
* `(upper | guards) - sum` similarly tests `sum <= upper`.

Only the guard bits corresponding to active `dl`/`du` constraints are checked. This is exact broadword arithmetic, not a hash or approximate filter.

Before high-level integerized mining, column minima are subtracted and the target is adjusted by `len*minimum`, exactly preserving every fixed-cardinality subset sum. All safety-sensitive multiplication/subtraction is checked; if an int64-safe transformation or packing plan cannot be created, the code falls back to the v0.2 engine.

`search_md_i8_packed` still uses the v0.2 completion-envelope values, but packs those envelope vectors once. The DFS loop then updates and checks `packed_lanes` integers rather than every original dimension.

### v0.3 PAT-inspired engine

`search_md_i8_pat` is used only for comonotonic matrices. A state stores lower and upper admissible index vectors. Since every column is nondecreasing, the sums at the lower and upper vectors are coordinatewise bounds on every combination represented by that box.

Before branching, the engine repeatedly:

1. raises each lower index by binary searching for the first value that can still meet all lower target constraints under a maximal completion;
2. lowers each upper index by binary searching for the last value that can still meet all upper constraints under a minimal completion;
3. propagates strict index ordering;
4. rejects an empty/impossible state.

It then branches at the position with the smallest nonzero `UB-LB` gap, closely following the important control idea in FLSSS `mPAT::grow()` while using ordinary Fortran arrays rather than the original pointer-bedded cache layout.

`engine='auto'` does not force PAT. It first checks whether a row ordering produces a comonotonic matrix and uses PAT only when a constrained target is within roughly 10% of a feasible edge and its tolerance is not excessively broad. Otherwise packed DFS is used.


### v0.4 mPAT engine and triangular cache

`flsss_mpat` moves the interval-state implementation closer to the original `mPATclass.hpp`/`mvalFindBound.hpp` architecture. Each state carries `LB`, `UB`, packed `sumLB`, and packed `sumUB`. Lower-bound finding asks whether a proposed index can still meet all lower constraints under the maximal completion represented by the current upper bounds; upper-bound finding uses the symmetric minimal completion. The two sweeps repeat until no bound changes. The next branch is the position with the smallest nonzero gap.

The original C++ engine accelerates these operations with triangular matrices `M[k][start]` containing sums of consecutive rows. v0.4 adds the same mathematical cache. For a candidate at position `p`, strict ordering implies that only a consecutive suffix of the affected upper prefix (or lower suffix) must change. The unchanged old contribution is obtained from a per-state packed prefix sum and the new consecutive contribution is fetched from the triangular cache. This changes each candidate sum construction from roughly `O(len*packed_lanes)` to `O(packed_lanes)` after the bound-prefix setup.

The Fortran cache is lane-major, `value(:,entry)`, which makes all lanes for one lookup contiguous. It contains `(2*N-len+1)*len/2` entries. Child states inherit their packed bound sums and update them incrementally during a split.

This is deliberately close to the computational ideas in `mPAT::grow()` and `findBoundCpp`, but it does not copy the original byte-bedded pointer layout or its residual/J-pointer search line-for-line. `engine='pat'` preserves the v0.3 implementation for direct A/B comparison.

`engine='auto'` now selects mPAT for moderate/large comonotonic searches, including central targets, while small or non-comonotonic problems continue to use packed DFS.

## Decomposition and OpenMP

`decompose_mflsss` and `decompose_arb_flsss` partition the admissible first-index range into exact, nonoverlapping contiguous ranges.

v0.3 adds `mflsss_decomp_run` and `arb_flsss_decomp_run`. Each decomposition object writes to its own result slot, so the OpenMP `parallel do schedule(dynamic)` loops require no locks during search. Results are merged deterministically by partition index afterward.

OpenMP is optional. `flsss_parallel::openmp_enabled()` uses the standard OpenMP conditional sentinel to distinguish an OpenMP build without creating a hard dependency on the OpenMP runtime in serial builds. Compatibility `maxCore` arguments activate the decomposed path only when OpenMP is present.

Serial decomposition shares one wall-clock budget across successive objects. In a parallel run the independent objects start under the same requested wall limit.

## Arbitrary precision and k-sum

Decimal strings are scaled columnwise to canonical signed arbitrary-length base-10 integers and are never converted to IEEE floating point.

The v0.2 k-sum accelerator remains unchanged in v0.3. Every k-combination stores exact indices, its exact sum tuple, a portable 64-bit fingerprint, and an ordering sorted by fingerprint. A matching fingerprint is always followed by full canonical arbitrary-precision tuple comparison, so collisions affect only performance.

## Knapsack and GAP

The exact knapsack and GAP reference implementations are unchanged in v0.4. They preserve the exported problems but do not attempt to reproduce every specialized C++ parallel upper-bound decomposition.

## Remaining differences from the original native engine

The largest remaining performance differences are:

* the original byte-bedded `PAT`/`mPAT` pointer/container layout and exact residual/J-pointer implementation inside `findBoundCpp`;
* original integer crunching/mask layouts beyond the portable guard-bit broadword representation;
* architecture-specific SIMD intrinsics and cache prefetch tuning;
* RcppParallel/TBB work stealing and distributed execution;
* original Bloom/xxHash organization (the Fortran arbitrary-precision lookup uses a sorted portable fingerprint index with exact verification);
* specialized GAP/knapsack threading and upper-bound decompositions.

## Fortran-specific safety

No algorithm relies on short-circuit evaluation of `.and.` or `.or.` when an index could be invalid. Recursive temporary state is invocation-local. Source lines stay within the standard 132-column free-form limit. All procedures are module procedures or have explicit interfaces.
