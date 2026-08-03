# Translation coverage

The upstream namespace exports 26 names. All are represented:

| R export | Fortran name | Status |
|---|---|---|
| `nloptr` | `nloptr` | Implemented |
| `nl.opts` | `nl_opts` | Implemented |
| `nloptr.get.default.options` | `nloptr_get_default_options` | Implemented |
| `nloptr.print.options` | `nloptr_print_options` | Implemented |
| `is.nloptr` | `is_nloptr` | Implemented |
| `nl.grad` | `nl_grad` | Implemented |
| `nl.jacobian` | `nl_jacobian` | Implemented |
| `check.derivatives` | `check_derivatives` | Implemented |
| `lbfgs` | `lbfgs` | Adapted backend |
| `varmetric` | `varmetric` | Adapted backend |
| `tnewton` | `tnewton` | Adapted backend |
| `neldermead` | `neldermead` | Implemented family |
| `sbplx` | `sbplx` | Adapted backend |
| `cobyla` | `cobyla` | Adapted backend |
| `bobyqa` | `bobyqa` | Adapted backend |
| `newuoa` | `newuoa` | Adapted backend |
| `slsqp` | `slsqp` | Adapted backend |
| `mma` | `mma` | Adapted backend |
| `ccsaq` | `ccsaq` | Adapted backend |
| `auglag` | `auglag` | Adapted backend |
| `direct` | `direct` | Adapted backend |
| `directL` | `direct_l` | Adapted backend |
| `crs2lm` | `crs2lm` | Adapted backend |
| `isres` | `isres` | Adapted backend |
| `stogo` | `stogo` | Adapted backend |
| `mlsl` | `mlsl` | Adapted backend |

The native NLopt C/C++ implementation and its language-binding/build machinery
are not translated line by line. The retained original wrapper sources document
the upstream interface and licensing.
