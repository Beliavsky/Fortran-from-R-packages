# Translation coverage

Upstream namespace exports two computational routines; both are translated.

| R export | Fortran API | Coverage |
|---|---|---|
| `solve.QP` | `solve_qp` | Complete |
| `solve.QP.compact` | `solve_qp_compact` | Complete |

The private Goldfarb-Idnani kernels, positive-definite factorization/solve, and
factor-inverse routines were also translated and modernized.

R registration code, `.Fortran` marshaling, package namespace machinery, and
R list construction are not applicable to the FPM library.
