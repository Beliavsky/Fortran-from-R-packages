# Notice

This project is derived from the computational interface and mathematical model of ECOSolveR and ECOS.

- ECOSolveR 0.6.1: GPL >= 3.
- Referenced ECOS `r-patches` source: GPL v3-or-later notices.
- This Fortran translation: GPL-3.0-or-later.

## Sparse linear algebra

The sparse symbolic/numeric `LDL^T` algorithm in `src/ecos_sparse.f90` follows the published SuiteSparse LDL algorithm used by ECOS.

SuiteSparse LDL is Copyright (c) Timothy A. Davis and is distributed under the GNU Lesser General Public License version 2.1 or, at the user's option, any later version. A copy of LGPL-2.1 is retained at:

```text
LICENSES/SuiteSparse-LDL-LGPL-2.1.txt
```

The Fortran module is distributed as part of this GPL-3.0-or-later combined work; the upstream LDL attribution and LGPL terms are retained here.

Version 0.4.0 also contains a native approximate-minimum-degree-style ordering implementation written for this translation. It is not a copy or line-for-line port of SuiteSparse AMD. RCM remains available for comparison/fallback.

## Version 0.4.0 certificate machinery

The sparse solver can construct post-failure homogeneous primal or dual ray problems to diagnose primal infeasibility or primal unboundedness without converting the original sparse problem to dense storage. This is not presented as the exact integrated `tau/kappa` homogeneous-self-dual embedding used by upstream ECOS; the distinction is documented in `TRANSLATION_COVERAGE.md`.

## Optional MatrixExtra adapter

`integration/matrixextra-adapter/vendor/MatrixExtra-fortran` is the previously translated MatrixExtra package and retains its own GPL-3.0-only notices and the licenses of its Matrix dependency. It is optional and is not linked into the core ECOS package. The adapter is intended for sparse construction/conversion interoperability; ECOS's KKT factorization remains in the core sparse backend.

The supplied ECOSolveR source archive is retained under `original/`.
