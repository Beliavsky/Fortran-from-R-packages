# API map

All 32 non-graphics NAMESPACE exports have Fortran counterparts.

| R export | Fortran mapping |
|---|---|
| `displ` | `displ` constructor |
| `conpl` | `conpl` constructor |
| `disexp` | `disexp` constructor |
| `conexp` | `conexp` constructor |
| `dislnorm` | `dislnorm` constructor |
| `conlnorm` | `conlnorm` constructor |
| `dispois` | `dispois` constructor |
| `conweibull` | `conweibull` constructor |
| `dist_pdf` | `dist_pdf`, `powerlaw_dist%pdf` |
| `dist_cdf` | `dist_cdf`, `powerlaw_dist%cdf` |
| `dist_all_cdf` | `dist_all_cdf` |
| `dist_data_cdf` | `dist_data_cdf` |
| `dist_data_all_cdf` | `dist_data_all_cdf` |
| `dist_ll` | `dist_ll`, `powerlaw_dist%loglik` |
| `dist_rand` | `dist_rand`, `powerlaw_dist%random` |
| `dpldis` | `dpldis` |
| `ppldis` | `ppldis` |
| `rpldis` | `rpldis` |
| `dplcon` | `dplcon` |
| `pplcon` | `pplcon` |
| `rplcon` | `rplcon` |
| `estimate_pars` | `estimate_pars` |
| `estimate_xmin` | `estimate_xmin` |
| `get_distance_statistic` | `get_distance_statistic` |
| `get_KS_statistic` | `get_KS_statistic` |
| `compare_distributions` | `compare_distributions` |
| `bootstrap` | `bootstrap` |
| `bootstrap_p` | `bootstrap_p` |
| `get_bootstrap_sims` | `get_bootstrap_sims` |
| `get_bootstrap_p_sims` | `get_bootstrap_p_sims` |
| `get_n` | `get_n`, `powerlaw_dist%get_n` |
| `get_ntail` | `get_ntail`, `powerlaw_dist%get_ntail` |

The exported R class names map to the `powerlaw_dist` derived type and its
family constants.  R's `plot`, `lines`, `points`, `show`, and the S3 bootstrap
plot methods are presentation-only and are not translated.
