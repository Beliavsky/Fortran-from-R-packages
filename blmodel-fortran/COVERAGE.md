# Computational coverage

## Exported upstream routines

| R routine | Fortran routine | Status |
|---|---|---|
| `BL_post_distr` | `bl_post_distribution` | Complete |
| `equilibrium_mean` | `equilibrium_mean` | Complete |
| `observ_normal` | `observ_normal` | Complete |
| `observ_powerexp` | `observ_powerexp` | Complete |
| `observ_ts` | `observ_student_t`, alias `observ_ts` | Complete |

## Internal upstream routines

| R routine | Fortran routine | Status |
|---|---|---|
| `.diag_of` | `diag_of` | Complete |
| `.make_diag` | `make_diag` | Complete |
| `.discrete_variance` | `discrete_variance` | Complete |
| `.equilibrium_mean_elliptic` | `equilibrium_mean_elliptic` | Complete |
| `.post_distr_new` | `post_distribution` | Complete |

## R-specific infrastructure

No numerical routine is omitted. R lists, matrix dimension names, `match.arg`,
dynamic `match.fun`, roxygen documentation, and testthat wrappers are replaced
by Fortran derived types, explicit validation, and procedure callbacks.

The original package contains no plotting, compiled C/C++, or bundled numerical
data files.
