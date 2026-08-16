# R-to-Fortran API map

The R package exports 51 names. `%>%` is an imported presentation/convenience
operator and is not applicable to Fortran. The remaining 50 computational
exports are represented.

| R export | Fortran | Notes |
|---|---|---|
| `dgkw`, `pgkw`, `qgkw`, `rgkw` | same names | GKw d/p/q/r |
| `llgkw`, `grgkw`, `hsgkw` | same names | negative log-likelihood, gradient, Hessian |
| `dbkw`, `pbkw`, `qbkw`, `rbkw` | same names | Beta-Kumaraswamy |
| `llbkw`, `grbkw`, `hsbkw` | same names | BKw likelihood derivatives |
| `dkkw`, `pkkw`, `qkkw`, `rkkw` | same names | Kumaraswamy-Kumaraswamy |
| `llkkw`, `grkkw`, `hskkw` | same names | KKw likelihood derivatives |
| `dekw`, `pekw`, `qekw`, `rekw` | same names | exponentiated Kumaraswamy |
| `llekw`, `grekw`, `hsekw` | same names | EKw likelihood derivatives |
| `dmc`, `pmc`, `qmc`, `rmc` | same names | McDonald family |
| `llmc`, `grmc`, `hsmc` | same names | McDonald likelihood derivatives |
| `dkw`, `pkw`, `qkw`, `rkw` | same names | Kumaraswamy |
| `llkw`, `grkw`, `hskw` | same names | Kumaraswamy likelihood derivatives |
| `dbeta_`, `pbeta_`, `qbeta_`, `rbeta_` | same names | upstream parameterization Beta(gamma, delta+1) |
| `llbeta`, `grbeta`, `hsbeta` | same names | beta likelihood derivatives |
| `gkwgetstartvalues` | same name | allocatable-vector Fortran function |

The density/CDF/quantile routines are elemental, so a scalar parameter can be
combined naturally with a conformable array argument. R's recycling of arrays
with different non-unit lengths is intentionally not reproduced; Fortran's
normal conformability rules are used instead.
