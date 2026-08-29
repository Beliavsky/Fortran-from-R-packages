# Notices and attribution

This distribution contains two licensing layers.

1. `src/lbfgsb3_core.f90` is a free-form module translation of the L-BFGS-B
   3.0 numerical code by Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis
   Morales, and collaborators. The supplied R package states that this code is
   available under the BSD 3-clause license. The exact license material shipped
   with the source package is retained in
   `LICENSES/LBFGSB-BSD-3-Clause.txt`.
2. `src/lbfgsb3.f90`, the FPM packaging, tests, examples, and documentation form
   the translated wrapper layer and are distributed under GPL-2.0-only,
   consistently with the `lbfgsb3c` package.

The original package identifies Matthew L. Fidler and John C. Nash as wrapper
and package authors. Original source files needed for provenance and comparison
are retained under `original_source/lbfgsb3c-main/`.

Publications describing work using L-BFGS-B should cite at least one of:

- R. H. Byrd, P. Lu, and J. Nocedal, "A Limited Memory Algorithm for Bound
  Constrained Optimization," SIAM Journal on Scientific Computing, 1995.
- C. Zhu, R. H. Byrd, and J. Nocedal, "L-BFGS-B: Algorithm 778," ACM
  Transactions on Mathematical Software, 1997.
- J. L. Morales and J. Nocedal, "Remark on Algorithm 778: L-BFGS-B," ACM
  Transactions on Mathematical Software, 2011.
