# API coverage notes

The authoritative function-by-function mapping is in `fpm.toml` under `[extra.translation]` and is mirrored in the README table.

## Storage types

- `csr_matrix`: values, one-based column indices, and one-based row pointers.
- `csc_matrix`: values, one-based row indices, and one-based column pointers.
- `coo_matrix`: values plus one-based row/column coordinate arrays.
- `ssr_matrix` / `ssc_matrix`: lower-triangle symmetric sparse storage.
- `csr_cholesky`: translated factorization object storing sparse CSR lower/upper factors plus Ng-Peyton ordering, elimination-tree, supernode, compressed `LINDX`/`XLINDX`, numerical `LNZ`/`XLNZ`, cache split, numerical level/tuning controls, tiny-pivot status, workspace, matrix determinant metadata, and SparseM-compatible triangular-factor determinant metadata.
- `csr_cholesky_stats`: `LSTATS`-style factor statistics returned by `chol_statistics`, including maximum supernode sizes and factor/solve operation counts.
- `slm_result` and `slm_summary_result`: numerical single-response linear-model results without R formula/class metadata.
- `mslm_result` and `mslm_summary_result`: multi-response fits and per-response numerical summaries.
- `hb_matrix`: real assembled Harwell-Boeing compressed-column data plus optional RHS, starting-guess, and exact-solution arrays.

## Remaining high-value parity work

1. Continue Ng-Peyton performance tuning on larger irregular supernodal problems. Cache-block partitioning and temporary-workspace sizing now match the upstream `FNSPLT`/`FNTSIZ` routines in direct oracle tests, and complete driver comparisons match permutation, symbolic storage, cache split, and `LNZ` exactly at levels 1/4/8 on tridiagonal, star, grid, and irregular sparse SPD cases. The modern level-2 kernel is kept mathematically consistent instead of reproducing divergent legacy behavior observed in the supplied source.
2. Expand `summary.slm`/`summary.mslm` parity with weighted-model metadata, residual summaries, and richer formula/call metadata where a non-R API is useful.
3. Expand R indexing/recycling compatibility for negative/logical subscripts and scalar/vector Ops.
4. Reduce dense temporary use in elementwise, binding, subsetting, and conversion paths.
