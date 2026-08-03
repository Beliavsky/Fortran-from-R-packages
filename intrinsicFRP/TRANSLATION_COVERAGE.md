# Translation coverage

| R export | Fortran procedure | Coverage |
|---|---|---|
| `TFRP` | `tfrp` | Implemented, including HAC errors |
| `FRP` | `frp` | FM and misspecification-robust variants |
| `SDFCoefficients` | `sdf_coefficients` | FM and GKR variants |
| `GKRFactorScreening` | `gkr_factor_screening` | Implemented |
| `OracleTFRP` | `oracle_tfrp` | GCV, CV, rolling, relaxed, errors |
| `FGXFactorsTest` | `fgx_factors_test` | Implemented with internal Lasso |
| `HJMisspecificationDistance` | `hj_misspecification_distance` | Point and interval estimates |
| `HACcovariance` | `hac_covariance` | Implemented with optional prewhitening |
| `IterativeKleibergenPaap2006BetaRankTest` | `iterative_kleibergen_paap_2006_beta_rank_test` | Adapted singular-value test |
| `ChenFang2019BetaRankTest` | `chen_fang_2019_beta_rank_test` | Deterministic bootstrap implementation |
| `GiglioXiu2021RiskPremia` | `giglio_xiu_2021_risk_premia` | Implemented with PCA selectors |

All 11 exported computational procedures are represented. Internal moment
estimators, adaptive weights, PCA selectors, HAC helpers, FGX covariance
routines, probability functions, random-number generation, and dense linear
algebra are also included.

Not compiled: R/Rcpp registration, data-frame/list conversion, package datasets,
and documentation-only infrastructure. The original computational R/C++ files
are retained under `original/`.
