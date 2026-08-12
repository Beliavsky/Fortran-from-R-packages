# Translation notes

## Upstream

- R package: `adagio`
- translated version: 0.9.2
- date in upstream DESCRIPTION: 2023-10-26
- author: Hans W. Borchers
- upstream license: GPL (>= 3)
- `NeedsCompilation: no`: the active package is pure R

The original attached package tree is retained under `original/adagio-master/`.

## Computational coverage

The following exported computational functions are represented:

`fnRosenbrock`, `grRosenbrock`, `fnRastrigin`, `grRastrigin`, `fnNesterov`,
`grNesterov`, `fnNesterov1`, `fnNesterov2`, `fnHald`, `grHald`, `fnShor`,
`grShor`, `fnTrefethen`, `fnWagon`, `maxquad`, `count`, `occurs`, `maxsub`,
`maxsub2d`, `maxempty`, `assignment`, `subsetsum`, `sss_test`, `setcover`,
`knapsack`, `mknapsack`, `change_making`, `bpp_approx`, `neldermead`,
`neldermeadb`, `hookejeeves`, `simpleEA`, `simpleDE`, `pureCMAES`,
`transfinite`, `hamiltonian`, and `Historize`.

`fminviz` and `flineviz` are omitted because they are plotting/interactive
visualization routines rather than numerical algorithms.

Unexported helper closures such as `initHald`, `initMaxquad`, `.hjsearch`, and
`graph_vect2list` are incorporated directly into the corresponding Fortran
implementations.

## lpSolve-backed routines

The original R source calls `lpSolve` in four places.  Earlier provisional work
used standalone replacements for these algorithms, but the final v0.1.0 port
instead uses the supplied `lpSolve-fortran v0.1.0` package as an FPM path
dependency.

- `assignment`: calls `lp_assign` and reconstructs the row permutation.
- `change_making`: minimizes a vector of ones subject to one equality, with all
  variables integer and nonnegative.
- `setcover`: minimizes set weights subject to element coverage constraints,
  with all decision variables binary.
- `mknapsack`: reproduces the R package's block binary MILP formulation,
  including the single-knapsack special case.

The vendored dependency is unmodified except for its location in the source
tree and retains its LGPL-2.0-only headers/license/provenance.

## Fortran representation choices

- R lists become explicit result derived types.
- R function closures become procedure callbacks or type-bound procedures.
- `Historize()` becomes `history_buffer`, whose `input(call,variable)` matrix
  stores one input vector per row.
- `transfinite()` returns two R closures; Fortran exposes those as
  `transfinite_forward` and `transfinite_inverse`.
- `maxquad()` returns two R closures over generated matrices; Fortran returns a
  `maxquad_problem` with type-bound `value` and `gradient` methods.
- Random algorithms use a standalone xorshift-style RNG and Box-Muller normal
  generation.  The same integer seed is reproducible within this Fortran port
  but is not expected to reproduce R's random stream.
- `pureCMAES` uses a self-contained symmetric Jacobi eigensolver rather than
  R's `eigen()`, keeping the package free of a BLAS/LAPACK requirement beyond
  the independent lpSolve dependency.

## Preserved source behavior

Several unusual source behaviors materially affect outputs and are kept:

1. `simpleDE` sets `nfeval <- N` initially, then adds `N` after each single
   objective evaluation.  `nfeval` therefore overcounts by a factor of `N` in
   the inner loop.  The Fortran result provides both `nfeval` and
   `actual_nfeval`.
2. `simpleEA` creates random initial candidates from `[0, upper-lower]`, without
   adding `lower`, exactly as written in the R code.  The same applies to its
   independent new candidates.
3. `simpleEA` evaluates the same source-defined row range
   `(N+1):(M*N+K)`, rather than all rows created by the `rbind` expression.
4. Hooke-Jeeves uses the source's `abs(fx) < target` continuation test.
5. Exact comparisons in `count` and `occurs` remain exact; they are not replaced
   by tolerance comparisons.

## Translation fixes required by Fortran semantics

R short-circuits `||` and `&&`; Fortran `.or.`/`.and.` do not guarantee
short-circuit evaluation.  Two literal translations were therefore unsafe and
were rewritten without changing the algorithm:

- unbounded Hooke-Jeeves no longer references unallocated `lower/upper` arrays;
- best/worst-fit bin packing no longer evaluates `b(0)` before the first
  candidate bin has been selected.

These are language-portability fixes, not changes to the R algorithms.

## Invalid-input handling

The R package frequently uses `stop()` for malformed inputs.  The Fortran API
cannot reproduce R exceptions directly; routines generally return empty result
arrays/default status fields for structurally invalid inputs.  Valid-input
numerical behavior is the compatibility target of v0.1.0.

## Small numerical adaptations

- Nelder-Mead uses `ceiling(maxfeval/100)` for the integer convergence-check
  interval, which reproduces the effective number of decrements of R's
  possibly fractional `kcount/100` value.
- `pure_cmaes` floors tiny/non-positive covariance eigenvalues at the positive
  floating-point minimum before forming inverse square roots.  This is a
  narrow roundoff safeguard around the R `sqrt(eigen(C)$values)` step.
