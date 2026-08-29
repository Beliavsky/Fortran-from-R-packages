# API mapping

This file maps the principal computational API of R `urca` 1.3-4 to the
Fortran routines in this port.

| R API | Fortran API | Module | Notes |
|---|---|---|---|
| `ur.df` | `adf_test` | `urca_unitroot` | Fixed/AIC/BIC lags; none/drift/trend |
| `ur.ers` | `ers_test` | `urca_unitroot` | DF-GLS and P-test; constant/trend |
| `ur.kpss` | `kpss_test` | `urca_unitroot` | Mu/tau stationarity forms |
| `ur.pp` | `pp_test` | `urca_unitroot` | Z-alpha/Z-tau and auxiliary statistics |
| `ur.sp` | `schmidt_phillips_test` | `urca_unitroot` | Tau/rho, polynomial degree 1:4 |
| `ur.za` | `zivot_andrews_test` | `urca_unitroot` | Intercept/trend/both break search |
| `ca.po` | `phillips_ouliaris` | `urca_cointegration` | Pu/Pz |
| `ca.jo` | `johansen_test` | `urca_cointegration` | Trace/eigen, deterministic/spec choices |
| `cajolst` | `johansen_level_shift` | `urca_breaks` | Endogenous level-shift search |
| `cajools` | `cajools_fit` | `urca_cointegration` | OLS representation from Johansen object |
| `alphaols` | `alphaols_fit` | `urca_cointegration` | Loading-matrix regression |
| `cajorls` | `cajorls_fit` | `urca_cointegration` | Rank-restricted VECM OLS |
| `blrtest` | `beta_restriction_test` | `urca_restrictions` | Restrictions on beta |
| `alrtest` | `alpha_restriction_test` | `urca_restrictions` | Restrictions on alpha |
| `ablrtest` | `alpha_beta_restriction_test` | `urca_restrictions` | Joint alpha/beta restrictions |
| `bh5lrtest` | `partly_known_beta_test` | `urca_restrictions` | Partly known cointegration space |
| `bh6lrtest` | `iterated_partly_known_beta_test` | `urca_restrictions` | Iterative partly known beta procedure |
| `lttest` | `linear_trend_lr_test` | `urca_restrictions` | Linear-trend restriction LR test |
| `punitroot` | `mackinnon_pvalue` | `urca_mackinnon` | MacKinnon response surfaces |
| `qunitroot` | `mackinnon_quantile` | `urca_mackinnon` | MacKinnon response surfaces |
| `unitrootTable` | `unitroot_table` | `urca_mackinnon` | Critical-value table generation |

## Important constants

The numerical APIs use integer named constants instead of character strings.
Import the module constants rather than hard-coding their values.

### ADF

- `UR_NONE`
- `UR_DRIFT`
- `UR_TREND`
- `LAG_FIXED`
- `LAG_AIC`
- `LAG_BIC`

### ERS

- `ERS_DFGLS`
- `ERS_PTEST`
- constant/trend model constants are exported by `urca_unitroot`

### KPSS

- `KPSS_MU`
- `KPSS_TAU`

### Phillips-Perron

- `PP_ZALPHA`
- `PP_ZTAU`
- `PP_CONSTANT`
- `PP_TREND`

### Schmidt-Phillips

- `SP_RHO`
- `SP_TAU`

### Zivot-Andrews

- `ZA_INTERCEPT`
- `ZA_TREND`
- `ZA_BOTH`

### Johansen

- `JO_EIGEN`
- `JO_TRACE`
- `JO_NONE`
- `JO_CONST`
- `JO_TREND`
- `JO_LONGRUN`
- `JO_TRANSITORY`

### Phillips-Ouliaris

- `PO_NONE`
- `PO_CONST`
- `PO_TREND`
- `PO_PU`
- `PO_PZ`

### MacKinnon

- `MACK_TAU`
- `MACK_NORM`
- `MACK_NC`
- `MACK_C`
- `MACK_CT`
- `MACK_CTT`

## Result objects

`urca_types` provides allocatable derived types rather than R S4 objects:

- `ur_test_result`
- `po_result`
- `johansen_result`
- `restriction_result`
- `vecm_result`
- `lm_result`

These carry numerical results such as statistics, critical values,
coefficients, residuals, eigenvalues/eigenvectors, covariance matrices,
break points, selected lags, p-values, and status codes.

The umbrella module `urca` re-exports the public modules so most clients can
simply write:

```fortran
use urca
```
