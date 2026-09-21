# funData - modern Fortran translation

This directory translates the computational core of the R package **funData** 1.3-9 into modern free-form Fortran with FPM packaging. The upstream package is by Clara Happ-Kurz and is distributed under GPL-2. This translation preserves that license and records upstream provenance in `NOTICE.md` and `upstream/`.

The Fortran API uses native derived types for regular, irregular, and multivariate functional data. Regular data are stored in column-major order with an explicit dimension vector and support grids. IEEE NaNs represent missing regular-grid observations. Numerical integration supports one- through three-dimensional regular supports, matching the upstream computational limit.

No external numerical dependency is required. The target `Fortran-from-R-packages` repository was checked before implementation; this translation does not copy `abind`, `fda`, BLAS, LAPACK, or any shared package source.

## Build

With FPM installed:

```text
fpm build
fpm test
fpm run --example basic_usage
```

GNU Fortran 14.2 direct-build validation is documented in `docs/VALIDATION.md`.

## Main API

`funData_api` re-exports the single `dp` kind, functional-data types, constructors needed by Fortran callers, and the translated computational operations. Important groups include:

- conversion between regular, irregular, and one-component multivariate representations;
- extraction/subsetting and IEEE-NaN interpolation;
- quadrature, integration, L2 norms, scalar products, pointwise means, and orientation flipping;
- two- and three-way tensor products;
- eigenfunction/eigenvalue families and univariate/multivariate simulation;
- sparsification and Gaussian error addition.

R-specific S4 dispatch, plotting, printing, summaries, data-frame conversion, names/accessor plumbing, and conversions that merely delegate to the external `fda` package are intentionally outside the computational coverage basis.

## Translation coverage

Package status: **substantial**. Coverage is **20 of 20 (100.0%)** exported computational R functions. This fraction measures mapped computational functions, not complete R compatibility. R S4/formal dispatch, object metadata, warning text, non-standard evaluation, package-specific RNG streams, and external-package conversion behavior are not implied by the percentage.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `.intWeights` | `funData_numeric` | `int_weights` | complete |
| `addError` | `funData_simulation` | `add_error_fun_data`, `add_error_multi_fun_data` | substantial |
| `approxNA` | `funData_numeric` | `approx_na` | substantial |
| `as.funData` | `funData_core` | `as_fun_data` | substantial |
| `as.irregFunData` | `funData_core` | `as_irreg_fun_data` | substantial |
| `as.multiFunData` | `funData_core` | `as_multi_fun_data` | substantial |
| `eFun` | `funData_simulation` | `eigenfunctions` | substantial |
| `eVal` | `funData_simulation` | `eigenvalues` | complete |
| `extractObs` | `funData_core` | `extract_fun_data`, `extract_irreg_fun_data`, `extract_multi_fun_data` | substantial |
| `flipFuns` | `funData_numeric` | `flip_fun_data`, `flip_fun_irreg_data`, `flip_irreg_fun_data`, `flip_multi_fun_data` | substantial |
| `integrate` | `funData_numeric` | `integrate_fun_data`, `integrate_irreg_fun_data`, `integrate_multi_fun_data` | substantial |
| `meanFunction` | `funData_numeric` | `mean_fun_data`, `mean_irreg_fun_data`, `mean_multi_fun_data` | substantial |
| `scalarProduct` | `funData_numeric` | `scalar_product_fun_data`, `scalar_product_fun_irreg`, `scalar_product_irreg_fun`, `scalar_product_irreg_fun_data`, `scalar_product_multi_fun_data` | substantial |
| `simFunData` | `funData_simulation` | `sim_fun_data` | substantial |
| `simMultiFunData` | `funData_simulation` | `sim_multi_fun_data` | substantial |
| `sparsify` | `funData_simulation` | `sparsify_fun_data`, `sparsify_multi_fun_data` | substantial |
| `tensorProduct` | `funData_numeric` | `tensor_product2`, `tensor_product3` | complete |
| `norm` | `funData_numeric` | `norm_fun_data`, `norm_irreg_fun_data`, `norm_multi_fun_data` | substantial |
| `subset` | `funData_core` | `extract_fun_data`, `extract_irreg_fun_data`, `extract_multi_fun_data` | substantial |
| `[` | `funData_core` | `extract_fun_data`, `extract_irreg_fun_data`, `extract_multi_fun_data` | substantial |

See `docs/API_COVERAGE.md` for exclusions and compatibility notes.

## License and citation

The upstream package declares GPL-2. The translation is distributed under the same license; see `LICENSE` and `NOTICE.md`. For academic use of the original package, retain the upstream citation in `upstream/inst/CITATION`.
