# API map

This map relates the attached R package's numerical routines to the modern Fortran procedures. R-only object, plotting, and slider wrappers are excluded.

## Archimedean

| R routine | Fortran procedure | Status |
|---|---|---|
| `archmList` | integer family IDs `1:22` | represented |
| `archmParam` | `archm_default_alpha` | tested |
| `archmRange` | `archm_range` | tested |
| `archmCheck` | `archm_check` | tested |
| `Phi`, `.Phi`, `.Phi0` | `archm_phi`, `archm_phi0` | tested for all 22 |
| `.invPhi` | `archm_inv_phi` | tested for all 22 |
| derivative helpers | `archm_phi_derivative`, `archm_inv_phi_derivative` | tested through densities/dependence calculations |
| `Kfunc`, `.Kfunc`, `.invK` | `archm_k`, `archm_inv_k` | tested through simulation |
| `parchmCopula` | `archm_cdf` | tested for all 22 |
| `darchmCopula` | `archm_density` | tested for all 22 |
| `rarchmCopula` | `archm_rng` | tested for all 22 |
| `pgumbelCopula`, `dgumbelCopula`, `rgumbelCopula` | family 4 through the generic CDF/density/RNG | tested |
| `archmTau`, `archmRho` | `archm_tau`, `archm_rho` | tested |
| `archmTailCoeff` | `archm_tail_coeff` | tested |
| `archmCopulaFit` | `archm_fit` | fit path tested for all 22 |
| `archmCopulaSim` | `archm_rng` plus `archm_fit` | numerical components tested |

## Elliptical

| R routine | Fortran procedure | Status |
|---|---|---|
| `ellipticalList` | names `norm`, `cauchy`, `t`, `logistic`, `laplace`, `kotz`, `epower` | represented |
| `ellipticalParam`, `ellipticalRange`, `ellipticalCheck` | fitter defaults and explicit bounds | represented/tested |
| `gfunc` | `elliptical_generator` | tested through densities |
| `.delliptical` | `elliptical_marginal_pdf` | tested for all seven |
| `.pelliptical` | `elliptical_marginal_cdf` | tested for all seven |
| `.qelliptical` | `elliptical_marginal_quantile` | inversion-tested for all seven |
| `dellipticalCopula` | `elliptical_copula_density` | tested for all seven |
| `pellipticalCopula` | `elliptical_copula_cdf` | tested for all seven |
| `rellipticalCopula` | `elliptical_rng` | tested for all seven |
| dedicated Normal/Cauchy/t wrappers | generic procedures with corresponding family name | tested |
| `ellipticalTau`, `ellipticalRho` | `elliptical_tau`, `elliptical_rho` | tested |
| `ellipticalTailCoeff` | `elliptical_tail_coeff` | tested |
| `ellipticalCopulaFit` | `elliptical_fit` | fit path tested for all seven |
| `ellipticalCopulaSim` | `elliptical_rng` plus `elliptical_fit` | numerical components tested |

## Extreme value

| R routine | Fortran procedure | Status |
|---|---|---|
| `evList` | five official family names | represented |
| `evParam` | `ev_default_param` | tested |
| `evRange` | `ev_bounds` | tested through validation/fitting |
| `evCheck` | `ev_check` | tested |
| `Afunc` and derivatives | `ev_dependence`, `ev_dependence_derivative` | tested |
| `pevCopula` | `ev_cdf` | tested for all five |
| `devCopula` | `ev_density` | tested for all five |
| `revCopula` | `ev_rng` | tested for all five |
| `evTau`, `evRho` | `ev_tau`, `ev_rho` | tested |
| `evTailCoeff` | `ev_tail_coeff` | tested |
| `evCopulaFit` | `ev_fit` | fit path tested for all five |
| `evCopulaSim` | `ev_rng` plus `ev_fit` | numerical components tested |
| supplemental `gumbelII`, `pi`, `m` Pickands cases | `ev_dependence`/`ev_cdf` | evaluation tested; no fitting claim |

## Empirical and utilities

| R routine | Fortran procedure | Status |
|---|---|---|
| `pempiricalCopula` | `empirical_copula_cdf`, `empirical_copula_grid` | tested |
| `dempiricalCopula` | `empirical_density_grid` | mass-tested |
| `.Debye`, `.Debye1` | `debye_function` | reference-tested |
| `pfrechetCopula` | `frechet_copula_cdf` | tested for `m`, `pi`, `w`, and `psp` |
| `.pmoCopula` | `marshall_olkin_cdf` | corrected standard formula, tested |
| `.dmoCopula` | not exposed | source routine duplicates probability and ignores singular density |
| rank conversion | `pseudo_observations` | used/tested through empirical paths |
| concordance calculations | `kendall_tau_sample`, `spearman_rho_sample` | used by fitting/initialization |

## Excluded infrastructure

All slider, contour, perspective, heatmap, S4 class, show/print, R attribute, list/grid adapter, and `timeSeries` metadata routines are excluded.
