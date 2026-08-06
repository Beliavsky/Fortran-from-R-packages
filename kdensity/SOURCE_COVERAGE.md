# Source coverage

| Upstream file | Coverage |
|---|---|
| `R/kdensity.R` | Core estimator, normalization, support checks, bandwidth handling, evaluation |
| `R/kernels.R` | All built-in numerical kernels and aliases |
| `R/bandwidths.R` | JH, RHE, HS, UCV, nrd0/nrd and compatibility BCV/SJ |
| `R/starts.R` | Common built-in starts, aliases, support metadata, custom callback model |
| `R/helpers.R` | Automatic support/kernel logic and compatibility checks |
| `R/generics.R` | Numerical `coef`/`logLik` content retained in result fields; graphics/S3 omitted |

No plotting code is translated. R environment mutation, expression parsing,
namespace discovery, and S3 replacement methods have no direct numerical
Fortran equivalent.
