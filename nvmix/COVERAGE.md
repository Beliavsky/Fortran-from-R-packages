# Computational coverage

Upstream version: **nvmix 0.1-2**.

## Public numerical families

| Upstream family | Fortran implementation | Status |
|---|---|---|
| `dnvmix`, `pnvmix`, `rnvmix`, `qnvmix` | Typed `nvmix_model` plus compatibility wrappers | Implemented |
| `dgnvmix`, `pgnvmix`, `rgnvmix` | Group mappings in `nvmix_model` | Implemented |
| `dNorm`, `pNorm`, `rNorm`, `fitNorm` | Closed-form density; typed probability/simulation/fitting | Implemented |
| `rNorm_sumconstr` | Exact orthogonal projection onto a weighted-sum hyperplane | Implemented |
| `dStudent`, `pStudent`, `rStudent`, `fitStudent` | Closed-form density, probability engine, simulation and profile-EM fit | Implemented |
| `dgStudent`, `pgStudent`, `rgStudent` | Independent inverse-gamma group mixtures | Implemented |
| Student copula functions | Marginal transforms and joint-density ratios | Implemented |
| Grouped Student copula functions | Group-specific marginal transforms | Implemented |
| Generic NVM copula functions | Model-derived marginal CDFs and quantiles | Implemented |
| `dgammamix`, `pgammamix`, `qgammamix`, `rgammamix` | Squared-Mahalanobis gamma mixtures | Implemented |
| `VaR_nvmix`, `ES_nvmix` | Exact normal/Student formulas and generic simulation | Implemented |
| `corgnvmix` | Moment-based grouped-mixture correlation | Implemented |
| Kendall/Spearman dependence | Elliptical closed forms | Implemented |
| `lambda_gStudent` | Exact equal-df formula; grouped approximation | Partial; documented |
| `rskewt`, `dskewt` | Shifted inverse-gamma normal mixture | Implemented |
| skew-t CDF and quantiles | Deterministic mixing integration and bisection | Implemented |
| skew-t copula | Density and simulation | Implemented |
| `fitnvmix` | Constant, inverse-gamma, Pareto and gamma families | Implemented with native estimators |
| `fitStudentcopula` | Profile search with transformed sample correlation | Implemented |
| `fitgStudentcopula` | Marginal tail estimates plus grouped correlation | Implemented |
| `qqplot_maha` | Sorted observed and theoretical distances | Numerical data implemented |

## Internal numerical support

- Normal, Student, gamma, chi-square and F density/CDF/quantile routines.
- Inverse-gamma density and quantiles.
- Deterministic RNG with optional seeds.
- Halton point generation.
- Cholesky factorization and SPD solves.
- Matrix inversion, covariance/correlation conversion and Jacobi eigenanalysis.
- Positive-semidefinite covariance repair support.

## Deliberately excluded

- R S3 `print`, `summary` and plotting methods.
- Graphical QQ plots and goodness-of-fit annotations.
- R data-frame, matrix-class and method-dispatch infrastructure.
- R dynamic-library registration and `.Call` wrappers.
- Bundled `.RData` files as compiled inputs.
- Exact adapters to `qrng`, `Matrix`, `copula`, `pcaPP`, `ADGofTest`, `mnormt` and `pracma`.
- Arbitrary R distribution names supplied through `get_mix_`.

Callers can add additional mixing laws by extending `nvmix_mixing` or by
constructing samples externally and using the exposed numerical components.
