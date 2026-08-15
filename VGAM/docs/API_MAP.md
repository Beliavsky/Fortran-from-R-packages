# Computational API map

| VGAM concept / R family | Fortran API |
|---|---|
| link functions (`logitlink`, `probitlink`, `clogloglink`, etc.) | `vgam_links` |
| normal/gamma/beta/binomial/Poisson helpers | `vgam_distributions` |
| `dirichlet` | `dirichlet_pdf`, `random_dirichlet`, `fit_dirichlet_regression` |
| Dirichlet expected information | `dirichlet_eim_shape`, `dirichlet_eim_logshape` |
| Gumbel/Frechet/Rayleigh/Pareto and related helpers | `vgam_distributions` |
| GEV / generalized Pareto | `vgam_extremes` |
| Kumaraswamy / Laplace | `vgam_extremes` |
| actuarial distributions | `vgam_actuarial` |
| special functions (`lambertW`, zeta, expint, polygamma) | `vgam_special` |
| `gaussianff` | `fit_gaussian` |
| `poissonff` | `fit_poisson` |
| `binomialff` | `fit_binomial` |
| `gammaff` | `fit_gamma` |
| `inverse.gaussianff` | `fit_inverse_gaussian` |
| `multinomial` baseline-category model | `fit_multinomial` |
| cumulative ordinal model | `fit_ordinal` |
| beta regression family | `fit_beta_regression` |
| negative-binomial family | `fit_negative_binomial` |
| zero-inflated Poisson | `fit_zero_inflated_poisson` |
| zero-truncated Poisson | `fit_zero_truncated_poisson`, `dztpois_v` |
| hurdle Poisson | `fit_hurdle_poisson`, `dhurdlepois_v` |
| hurdle negative binomial | `fit_hurdle_negative_binomial`, `dhurdlenbinom_v` |
| zero-inflated negative binomial | `fit_zero_inflated_negative_binomial`, `dzinb_v` |
| VGLM coefficient constraints / parallel effects | `fit_constrained_vglm`, `parallel_constraint` |
| `rrvglm` core alternating reduced-rank model | `fit_rrvglm`, `rrvglm_result_t` |
| `drrvglm` H.A/H.C constrained reduced-rank model | `fit_drrvglm`, `drrvglm_result_t` |
| `qrrvglm` quadratic reduced-rank core | `fit_qrrvglm`, `qrrvglm_result_t` |
| `cqo` constrained quadratic ordination core | `fit_cqo`, `cqo_response_surface`, `cqo_calibrate` |
| `cao` rank-1 constrained additive ordination | `fit_cao_rank1`, `cao_result_t` |
| GAITD discrete alteration concepts | `gaitd_transform_pmf`, `gaitd_poisson`, `gaitd_negative_binomial` |
| GAITD Poisson/NB special-mass regression | `fit_gaitd_poisson_regression`, `fit_gaitd_nb_regression` |
| GAITD MLM direct alteration/inflation/deflation | `gaitd_mlm_poisson`, `gaitd_mlm_negative_binomial` |
| GAITD MLM regression | `fit_gaitd_mlm_poisson_regression`, `fit_gaitd_mlm_nb_regression` |
| GAITD outer-distribution mix semantics | `gaitd_mix_poisson`, `gaitd_mix_negative_binomial` |
| GAITD outer-mix regression | `fit_gaitd_mix_poisson_regression`, `fit_gaitd_mix_nb_regression` |
| GAITD NB outer-mix separate dispersion regression | `fit_gaitd_mix_nb_dispersion_regression` |
| Clayton copula | `clayton_copula_pdf`, `clayton_copula_cdf`, `random_clayton_copula` |
| Frank copula | `frank_copula_pdf`, `frank_copula_cdf`, `random_frank_copula` |
| FGM copula | `fgm_copula_pdf`, `fgm_copula_cdf`, `random_fgm_copula` |
| Gaussian copula | `gaussian_copula_pdf`, `gaussian_copula_cdf`, `random_gaussian_copula` |
| Plackett copula | `plackett_copula_pdf`, `plackett_copula_cdf`, `random_plackett_copula` |
| Ali-Mikhail-Haq copula | `amh_copula_pdf`, `amh_copula_cdf`, `random_amh_copula` |
| univariate Student-t helpers | `student_t_pdf`, `student_t_cdf`, `student_t_quantile`, `random_student_t` |
| `bistudentt` / `dbistudentt` | `bivariate_student_t_pdf`, `fit_bivariate_student_t` |
| `dbistudenttcop` | `student_t_copula_pdf`, `random_student_t_copula`, `fit_student_t_copula` |
| `binormal` / `dbinorm` | `bivariate_normal_pdf`, `random_bivariate_normal`, `fit_bivariate_normal` |
| `bilogistic` / `dbilogis` / `pbilogis` | `bivariate_logistic_pdf`, `bivariate_logistic_cdf`, `fit_bivariate_logistic` |
| `freund61` | `freund61_pdf`, `random_freund61`, `fit_freund61` |
| `trinormal` | `trivariate_normal_pdf`, `random_trivariate_normal`, `fit_trivariate_normal` |
| `bifgmexp` | `bifgm_exponential_pdf`, `bifgm_exponential_cdf`, `fit_bifgm_exponential` |
| `cens.normal` | `fit_censored_normal`, `censored_normal_logprob` |
| `cens.poisson` | `fit_censored_poisson`, `censored_poisson_logprob` |
| `cens.exponential` | `fit_censored_exponential`, `censored_exponential_logprob` |
| censored Rayleigh | `fit_censored_rayleigh`, `censored_rayleigh_logprob` |
| `tobit` / `d/p/q/rtobit` | `fit_tobit`, `dtobit_v`, `ptobit_v`, `qtobit_v`, `rtobit_v` |
| `foldnormal` / `d/p/q/rfoldnorm` | `fit_folded_normal`, `dfoldnorm_v`, `pfoldnorm_v`, `qfoldnorm_v`, `rfoldnorm_v` |
| EIM/constraint information helpers | `score_outer_information`, `observed_information`, `constrained_information` |
| `d/p/q/rzoabeta` | `dzoabeta`, `pzoabeta`, `qzoabeta`, `rzoabeta` |
| zero/one-altered beta regression | `fit_zoa_beta_regression`, `zoa_beta_result_t` |
| positive normal | `dposnorm_v`, `pposnorm_v`, `qposnorm_v`, `rposnorm_v` |
| positive geometric | `dposgeom_v`, `pposgeom_v`, `qposgeom_v`, `rposgeom_v` |
| positive/zero-truncated Poisson | `dpospois_v`, `ppospois_v`, `qpospois_v`, `rpospois_v` |
| positive negative binomial | `dposnbinom_v`, `pposnbinom_v`, `qposnbinom_v`, `rposnbinom_v`, `fit_positive_negative_binomial` |
| zero-altered Poisson | `dzapois_v`, `pzapois_v`, `qzapois_v`, `rzapois_v`, `fit_zero_altered_poisson` |
| zero-altered negative binomial | `dzanbinom_v`, `pzanbinom_v`, `qzanbinom_v`, `rzanbinom_v`, `fit_zero_altered_negative_binomial` |
| zero-altered geometric/binomial | `dzageom_v`, `dzabinom_v`, `fit_zero_altered_geometric`, `fit_zero_altered_binomial` |
| zero-inflated/deflated geometric/binomial | `dzigeom_v`, `pzigeom_v`, `dzibinom_v`, `pzibinom_v` |
| zero/one-inflated beta-binomial helpers | `dzoibetabinom_ab`, `pzoibetabinom_ab`, `rzoibetabinom_ab` |
| copula dependence regression | `fit_copula_regression`, `copula_regression_result_t` |
| `yeo.johnson` | `yeo_johnson`, `yeo_johnson_inverse`, `dyj_dy` |
| Yeo-Johnson normal/LMS-style likelihood | `fit_yj_normal` |
| `lms.yjn` three-predictor likelihood core | `fit_lms_yj`, `lms_yj_result_t` |
| `dAR1` / Gaussian `AR1` core | `dar1_log`, `dar1_density`, `fit_ar1` |
| `garma` with upstream-supported `q.ma.lag=0` | `fit_garma`, `garma_result_t` |
| `rrar` nested reduced-rank autoregression | `fit_rrar`, `rrar_result_t` |
| `s()` / spline additive fitting concept | `fit_pspline_vglm`, `fit_gam_*` |
| B-spline / natural spline bases | vendored `splines` module |
| family optimizer helpers | `vgam_optim` |
| dense fitting algebra | `vgam_linalg` |

Use `use vgam` for the aggregate public API, or import individual modules to
keep namespaces small.
