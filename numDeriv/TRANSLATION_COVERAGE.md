# Translation coverage

| R export | Fortran API | Status |
|---|---|---|
| `grad` | `grad`, `grad_elementwise`, `grad_complex`, `grad_elementwise_complex` | Complete |
| `jacobian` | `jacobian`, `jacobian_complex` | Complete |
| `hessian` | `hessian`, `hessian_complex` | Complete |
| `genD` | `genD` / `gend` | Complete |

The split complex and elementwise entry points make the R routine's dynamic
result-type cases explicit and type-safe in Fortran.
