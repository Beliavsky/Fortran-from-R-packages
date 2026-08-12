# FLSSS-fortran 0.4.0

Modern Fortran/FPM translation of the computational surface exported by **FLSSS 9.2.8**.

FLSSS is a specialized exact-optimization package for subset-sum-family search, multidimensional knapsack and generalized assignment. The original implementation is highly engineered C++/Rcpp/RcppParallel code. This project preserves the exported mathematical problems in portable modern Fortran and progressively ports the important search accelerators while keeping simpler reference engines available for validation.

## Main API

Use:

```fortran
use flsss_mod
```

Preferred underscore spellings and R-name-compatible spellings are both available.

| R function | Preferred Fortran routine | Compatibility spelling |
|---|---|---|
| `FLSSS` | `flsss` | `flsss` |
| `FLSSSmultiset` | `flsss_multiset` | `flsssmultiset` |
| `mFLSSSpar` | `mflsss_par` | `mflssspar` |
| `mFLSSSparImposeBounds` | `mflsss_par_impose_bounds` | `mflsssparimposebounds` |
| `mFLSSSparIntegerized` | `mflsss_par_integerized` | `mflsssparintegerized` |
| `mFLSSSparImposeBoundsIntegerized` | `mflsss_par_impose_bounds_integerized` | `mflsssparimposeboundsintegerized` |
| `decomposeMflsss` | `decompose_mflsss` | `decomposemflsss` |
| `mFLSSSobjRun` | `mflsss_obj_run` | `mflsssobjrun` |
| `mmKnapsack` | `mm_knapsack` | `mmknapsack` |
| `mmKnapsackIntegerized` | `mm_knapsack_integerized` | `mmknapsackintegerized` |
| `GAP` | `gap_solve` | `gap` |
| `auxKnapsack01bb` | `aux_knapsack01bb` | `auxknapsack01bb` |
| `auxKnapsack01dp` | `aux_knapsack01dp` | `auxknapsack01dp` |
| `auxGAPbb` | `aux_gap_bb` | `auxgapbb` |
| `auxGAPbbDp` | `aux_gap_bbdp` | `auxgapbbdp` |
| `auxGAPga` | `aux_gap_ga` | `auxgapga` |
| `arbFLSSS` | `arb_flsss` | `arbflsss` |
| `decomposeArbFLSSS` | `decompose_arb_flsss` | `decomposearbflsss` |
| `arbFLSSSobjRun` | `arb_flsss_obj_run` | `arbflsssobjrun` |
| `ksumHash` | `build_ksum_hash` | `ksumhash` |
| `addNumStrings` | `add_num_strings` | `addnumstrings` |

Additional Fortran orchestration routines are `mflsss_decomp_run`, `arb_flsss_decomp_run`, and `mflsss_par_integerized_parallel`.

## Search capabilities

* fixed-size and variable-size 1-D subset sum;
* per-position lower/upper index bounds and conjugate search;
* multiset subset sum with requested cardinalities per bucket;
* multidimensional subset sum with FLSSS `dl`/`du` lower/upper-dimension semantics;
* decomposition into resumable exact subproblems;
* FLSSS-style integerization;
* exact arbitrary-precision decimal-string subset sum;
* active k-sum fingerprint lookup with exact arbitrary-precision collision verification;
* exact completion-envelope pruning for unsorted data;
* packed int64 multidimensional comparisons;
* native mPAT interval-state search with triangular packed-sum caching for comonotonic integerized data;
* optional OpenMP execution of decomposed searches.

## Integerized search engines

`mflsss_par_integerized` and `mflsss_par_impose_bounds_integerized` accept an optional `engine` argument:

```fortran
r = mflsss_par_integerized(len, v, target, me, engine='auto')
```

Supported values are:

* `"auto"` - default. Uses mPAT for moderate/large comonotonic searches, including central targets, and packed DFS otherwise.
* `"mpat"` - force the v0.4 triangular-cache mPAT engine; unsuitable row order safely falls back.
* `"pat"` - force the v0.3 interval-state reference engine.
* `"packed"` - v0.3 broadword packed DFS.
* `"dfs"` - v0.2 native-int64 completion-envelope engine, retained as a numerical reference.

The returned `subset_solutions` reports `engine`, `packed_lanes`, `nodes`, `pruned`, `bound_states`, `tri_entries`, `tri_lookups`, `bound_updates`, and `mpat_splits`.

### Broadword packing

For fixed-cardinality integerized problems, each dimension is first shifted to a nonnegative range by subtracting its column minimum and adjusting the target by exactly `len*minimum`. Several dimensions are then placed into one signed-64-bit lane. Every field has a guard bit and enough data bits for the largest possible subset sum, so row additions cannot carry into neighboring fields.

Lower and upper comparisons use guard-bit subtraction, allowing several independent dimension inequalities to be checked at once. This is exact: packing is used only when the required field widths fit safely; otherwise the code falls back to the v0.2 int64 engine.

### v0.4 mPAT and triangular sums

For comonotonic rows, `engine='mpat'` keeps the mPAT state variables that matter computationally: lower/upper admissible index vectors and their packed sums. Bound finding alternates lower and upper tightening to a fixed point and then bisects the coordinate with the smallest nonzero `UB-LB` gap, following the control structure of FLSSS `mPAT::grow()`/`update()`.

The v0.4 difference is a triangular consecutive-sum cache analogous to the original `TriM/M` structure. For every length `k=1..len` and feasible row start, it stores the packed sum of those `k` consecutive rows. When testing a candidate lower bound, only a consecutive suffix of the current maximal prefix changes; the old contribution comes from a cached per-state prefix sum and the replacement is one triangular lookup. Upper-bound tests use the symmetric construction. Child states carry their `sumLB`/`sumUB` values through the split instead of resumming selected rows.

The cache is stored lane-major in Fortran (`value(:,entry)`) so all packed lanes for one triangular entry are contiguous. The number of triangular entries is `(2*N-len+1)*len/2`. `engine='pat'` remains available as the v0.3 implementation for A/B comparison.

Unlike v0.3, `auto` is willing to use mPAT on central targets when the matrix is comonotonic and the problem is large enough for interval-state tightening to pay for its setup cost.

## Parallel decomposition

The source contains standard OpenMP sentinel directives. A normal build requires no OpenMP runtime and remains serial.

Serial build:

```text
fpm build
fpm test
```

OpenMP build with GNU Fortran:

```text
fpm build --flag "-fopenmp"
fpm test  --flag "-fopenmp"
```

Then `mflsss_decomp_run(..., parallel=.true., max_threads=4)` and `arb_flsss_decomp_run(...)` execute independent decomposition objects concurrently. R-name compatibility wrappers also honor `maxCore` when OpenMP is actually enabled; without `-fopenmp` they stay on the ordinary serial path.

## Example

```fortran
program demo
  use flsss_mod
  implicit none
  type(subset_solutions) :: ans
  real(dp) :: v(10)
  integer :: i

  v = [(real(i,dp), i=1,10)]
  ans = flsss(3, v, target=15.0_dp, me=0.0_dp, solution_need=5)

  do i=1,ans%size()
    print '(*(i0,1x))', ans%sol(i)%idx
  end do
end program demo
```

## v0.4 benchmark summary

GNU Fortran 14.2.0, `-O2`, local container measurements:

* Central 45-row, 24-dimensional, subset-size-10 exact search: packed DFS **~0.64 s / 4,168,901 states**, v0.3 PAT **~0.059 s / 20,274 states**, v0.4 mPAT **~0.045 s / 20,274 states**. This is about **14x** versus packed DFS and **1.3x** versus the v0.3 PAT engine at identical interval-state count. The mPAT run used 405 triangular entries and about 708k triangular lookups.
* The v0.3 broadword benchmark remains about **3x** faster than scalar int64 DFS when the search tree is unchanged.
* The v0.2 arbitrary-precision k-sum accelerator remains active and continues to give order-of-magnitude gains over plain arbitrary-precision DFS.
* OpenMP decomposition remains available unchanged from v0.3.

These are algorithm/workload demonstrations, not universal speed guarantees. See `VALIDATION.md` and `example/benchmark_mpat.f90`.

## Remaining performance differences from original FLSSS

v0.4 ports the main mPAT interval-state and triangular-sum ideas, but it is not a line-for-line copy of the native C++ engine. Remaining differences include the original byte-bedded pointer/container layout, the exact residual/J-pointer implementation inside `findBoundCpp`, architecture-specific SIMD and prefetch intrinsics, specialized multi-column mask layouts beyond the portable guard-bit representation, TBB work stealing/distributed execution, and xxHash/Bloom-specific hash-table organization.

The Fortran implementation remains exact on represented problems and retains the v0.2 engine for cross-checking.

## Licensing

FLSSS is GPL-3.0-only. See `LICENSES.md`, `LICENSE`, and `original/`.
