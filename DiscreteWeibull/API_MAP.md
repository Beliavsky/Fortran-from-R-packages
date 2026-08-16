# API map

The upstream package uses `exportPattern("^[[:alpha:]]+")`. The public
computational objects are mapped as follows.

| R API | Fortran API |
|---|---|
| `ddweibull` | `ddweibull` |
| `pdweibull` | `pdweibull` |
| `qdweibull` | `qdweibull` |
| `rdweibull` | `rdweibull` |
| `Edweibull` | `Edweibull` |
| `E2dweibull` | `E2dweibull` |
| `Vdweibull` | `Vdweibull` |
| `ERdweibull` | `ERdweibull` |
| `loglikedw` | `loglikedw` |
| `lossdw` | `lossdw` |
| `estdweibull` | `estdweibull` returning `dweibull_fit_result` |
| `varFisher` | `varFisher` returning `fisher_result` |
| `ddweibull3` | `ddweibull3` |
| `pdweibull3` | `pdweibull3` |
| `qdweibull3` | `qdweibull3` |
| `rdweibull3` | `rdweibull3` |
| `hdweibull3` | `hdweibull3` |
| `Edweibull3` | `Edweibull3` |
| `E2dweibull3` | `E2dweibull3` |
| `loglikedw3` | `loglikedw3` |
| `lossdw3` | `lossdw3` |
| `estdweibull3` | `estdweibull3` returning `dweibull_fit_result` |

The `Dq`, `Dbeta`, `Dqbeta`, `Dqq`, and `Dbetabeta` names appearing in
`varFisher.R` are local nested helper closures and are not package-level
exports. The Fortran `varFisher` computes the same observed information from
the fitted log-likelihood Hessian.
