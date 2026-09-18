# API coverage

Coverage is measured using distinct **exported computational R functions** in
`NAMESPACE`, following the repository translation convention. R6 methods that
form the behavior of the exported `varmapack_model()` constructor are described
under that one exported function rather than counted as additional R exports.

| R function | Fortran module | Fortran procedure(s) | Status |
|---|---|---|---|
| `varmapack_autocov` | `varmapack_analysis` | `varmapack_autocov` | complete |
| `varmapack_cov2corr` | `varmapack_analysis` | `varmapack_cov2corr` | complete |
| `varmapack_model` | `varmapack_model_mod` | `make_varmapack_model`, `model_sim`, `model_acvf`, `model_psi`, `model_irf`, `model_specrad`, `model_ma_specrad` | substantial |
| `varmapack_testcase` | `varmapack_testcases_mod` | `testcase_by_name`, `testcase_by_index` | substantial |
| `varmapack_testcases` | `varmapack_testcases_mod` | `varmapack_testcases` | substantial |

The model mapping includes:

- Gaussian VAR/VMA/VARMA simulation without burn-in for stationary models;
- conditional simulation from supplied startup observations;
- nonstationary simulation when startup values are supplied;
- VARMAX simulation with shared or replicate-specific exogenous paths;
- time-varying mean paths for VARMA;
- theoretical autocovariances;
- PSI and orthogonalized impulse responses;
- AR and MA spectral radii;
- all upstream named and parameterized testcase constructors.

See `PROVENANCE.md` for numerical and interface differences that prevent this
100% exported-function mapping fraction from implying exact R compatibility.
