# Rcsdp-fortran 0.3.0

Modern Fortran/FPM translation of the computational code in **Rcsdp 0.1.57.6**, including the bundled **CSDP 6.1.1** semidefinite-programming solver.

The R `.Call` interface and R class/coercion machinery are intentionally omitted. The numerical solver, sparse SDP data structures, SDPA I/O, and CSDP computational algorithms are provided directly as Fortran modules. The original Common Public License 1.0 and provenance files are preserved.

## v0.3.0 highlights

v0.3.0 extends the v0.2 sparse Schur implementation into the predictor/corrector block products and line search.

New in this release:

- native translation of CSDP `makefill.c`;
- native translations of the three restricted products from `mat_multsp.c`:
  - `mat_multspa`: A is restricted to fill,
  - `mat_multspb`: B is restricted to fill,
  - `mat_multspc`: only fill entries of C are generated;
- row- and column-indexed fill structures for efficient traversal;
- race-free OpenMP directives on restricted-product loops; the same sources compile serially without OpenMP;
- predictor reuse of `Zi*Fd`, following the original CSDP dataflow and avoiding redundant dense products;
- fill-restricted corrector products and true-operator Schur refinement;
- CSDP-style Schur diagonal equilibration (`workvec8` concept) and retry regularization;
- large-block Lanczos line search with full double reorthogonalization, while retaining exact LAPACK eigenvalue search for smaller blocks;
- retained v0.2/v0.1 fallback controls for regression testing;
- fill-product, solver-mode, and Lanczos regression tests;
- a dedicated fill-product benchmark.

## Solver functionality

The port supports:

- dense symmetric PSD blocks and diagonal/nonnegative blocks;
- sparse symmetric constraint blocks;
- sparse `A(X)` and `A'(y)` operators;
- CSDP primal-dual predictor-corrector Newton equations;
- sparse-aware Schur assembly

  `O(d) = A(Z^{-1} A'(d) X)`;

- dense LAPACK Cholesky of the Schur matrix;
- CSDP-style Schur equilibration and diagonal retry regularization;
- optional true-operator iterative refinement;
- fill-restricted predictor/corrector products;
- exact or Lanczos positive-definite line search;
- CSDP initialization, stopping measures, status codes, and controls;
- SDPA sparse problem and solution I/O;
- Rcsdp-style symmetric triplet helpers.

The public facade is `use rcsdp`.

## New controls in v0.3.0

`type(csdp_control)` includes:

```fortran
logical  :: use_sparse_schur       = .true.
logical  :: use_fill_products      = .true.
real(dp) :: fill_density_limit     = 0.01_dp
logical  :: use_schur_scaling      = .true.
logical  :: use_lanczos_linesearch = .true.
integer  :: lanczos_threshold      = 180
integer  :: lanczos_iterations     = 30
logical  :: fastmode               = .false.
```

`use_fill_products=.false.` retains the v0.2 block-product path. `use_sparse_schur=.false.` retains the v0.1 dense-reference Schur path. This makes the optimizations independently testable.

The 1% default fill-density cutoff follows the intent of CSDP's `SPARSELIMA/B/C` defaults. Blocks above the cutoff automatically use ordinary dense `matmul`.

For matrix blocks larger than `lanczos_threshold`, the line search estimates the limiting eigenvalue with at most `lanczos_iterations` Lanczos iterations. Smaller blocks use the exact LAPACK symmetric eigensolver.

## Solution diagnostics

In addition to the v0.2 Schur statistics, `csdp_solution` reports:

```text
fill_nnz
fill_full_entries
fill_sparse_products
fill_dense_products
lanczos_linesearches
schur_diagadd
```

These are useful for determining whether a particular SDP actually benefits from fill restriction.

## Build

Requirements: a Fortran 2018 compiler, BLAS, LAPACK, and optionally OpenMP.

```text
fpm build
fpm test
fpm run --example basic_example
fpm run --example solve_sdpa -- data/theta1.dat-s theta1.sol
fpm run --example benchmark_theta
fpm run --example benchmark_fill_products
```

GNU Fortran OpenMP build, if desired:

```text
fpm build --flag "-O3 -fopenmp"
fpm test  --flag "-O3 -fopenmp"
```

The `!$omp` directives are ignored by a normal serial build.

## Validation

The bounds-checked regression suite covers:

1. the canonical CSDP example (objective `2.75`);
2. diagonal-cone LP functionality;
3. Rcsdp's manifold example;
4. triplet and SDPA I/O round trips;
5. CSDP `theta1.dat-s` (objective `23.0`);
6. sparse-vs-dense `A`, `A'`, and Schur equality;
7. `spA`, `spB`, and `spC` against full dense products;
8. true-operator fill-restricted Schur application;
9. default v0.3, v0.2-style no-fill, and unscaled-Schur solver modes;
10. a 220x220 Lanczos line-search regression against the exact LAPACK result;
11. an integrated 120x120 sparse-fill SDP whose solver path exercises only restricted products and converges to objective `-1`.

The suite passes with GNU Fortran 14.2.0 under `-O0 -fcheck=all` and under `-O2`. The fill and `theta1` tests also pass with `-fopenmp` and two threads.

Illustrative local `-O2` results:

```text
theta1 sparse solver:       about 0.019-0.023 s
theta1 dense reference:     about 0.420-0.424 s
ratio:                       about 18-22x
```

For the included synthetic 500x500 fill benchmark (740 required entries out of 250,000, or 0.296% fill):

```text
full dense product:          about 0.0155-0.0157 s
fill-restricted spC:         about 0.00036-0.00037 s
ratio:                       about 42-44x
max required-entry error:    about 2-3e-13
```

Timings depend strongly on compiler, BLAS, CPU, threading, and problem sparsity; they are engineering regression numbers rather than performance guarantees.

## Remaining differences from CSDP 6.1.1

The main remaining low-level differences are narrower than in v0.2:

- CSDP packed temporary storage is replaced by ordinary allocatable Fortran arrays;
- compiler-specific prefetch builtins are not reproduced;
- OpenMP scheduling is simpler and portable rather than matching every original C pragma;
- the Lanczos routine uses a small LAPACK eigensolve of its tridiagonal projection instead of the original `qreig.c` Algorithm 464 implementation;
- some historical CSDP convergence/refinement heuristics are expressed more directly in modern Fortran;
- the Schur matrix remains dense, as in CSDP's direct solver architecture.

R S3/S4 dispatch, `Matrix` coercions, `.Call` registration, and R control-file glue remain intentionally omitted. There is no plotting code to translate.

## License and provenance

Rcsdp and bundled CSDP 6.1.1 declare **Common Public License Version 1.0 (CPL-1.0)**. The full license is retained as `LICENSE`; original metadata and authorship files are under `licenses/`.
