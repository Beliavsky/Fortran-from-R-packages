# Export coverage

The R namespace contains 279 export declarations. They are represented in the
Fortran project by category rather than by a one-to-one wrapper count.

| R export category | Fortran representation |
|---|---|
| ECD constructors, roots, density, CDF, moments, transforms | `ecd_core` |
| ECLD/SGED constructors, moments, transforms, option operators | `ecld_models` |
| Laplace, stable-count, SLD, QSLD, Levy distributions | `ecd_processes` |
| LAMP generation and simulation | `lamp_process` |
| Black-Scholes and option interpolation | `ecd_options` |
| Fitting and optimization | `ecd_fitting` |
| Lags, returns, quantiles, and tail statistics | `ecd_timeseries` |
| Special functions, quadrature, and roots | `ecd_math` |
| Short original-compatible names | `ecd_compat` |
| MPFR/vector/S4 dispatch wrappers | common typed procedures |
| Plot-only methods | numeric inputs already returned by core routines |
| Database/configuration/data loaders | omitted as non-numerical infrastructure |

The original package source is included so every R wrapper remains traceable.
