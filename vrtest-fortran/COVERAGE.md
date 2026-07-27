# Computational coverage

This table maps each R source file in vrtest 1.2 to its Fortran implementation.

| R source | Fortran implementation | Status |
|---|---|---|
| `ABEL1Q.R` | `abel_bandwidth` | Ported |
| `AR1.R` | `ar1_fit` | Ported; optional legacy indexing |
| `Adjust.thin.R` | `adjust_thin` | Ported |
| `Auto.Q.R` | `automatic_portmanteau` | Ported |
| `Auto.VR.R` | `automatic_variance_ratio` | Ported |
| `AutoBoot.test.R` | `automatic_vr_bootstrap` | Ported |
| `Ave.Ex.R` | `average_exponential_test` | Ported |
| `Boot.test.R` | `variance_ratio_bootstrap` | Ported |
| `Chen.Deo.R` | `chen_deo_test` | Ported |
| `Chow.Denning.R` | `chow_denning` | Ported |
| `DL.test.R`, `DLtest.R` | `dominguez_lobato_test`, `dominguez_lobato_statistic` | Ported |
| `FastLMCD.R`, `LMCD.R` | `lmcd_statistics` | Ported |
| `FastLM_stat.R`, `LM_stat.R` | `lm_statistic` | Ported |
| `FastVR.R` | `fast_variance_ratio` | Ported |
| `Gen.Spec.Test.R` | `generalized_spectral_test` | Ported |
| `IACM.R` | internal `integrated_acm` | Ported |
| `ISIMP.R` | `simpson_integral` | Ported |
| `JWright.crit.R` | `joint_wright_critical_values` | Ported |
| `Joint.Wright.R` | `joint_wright_tests` | Ported |
| `Lo.Mac.R` | `lo_mackinlay` | Ported |
| `Mammen.R`, `Rademacher.R` | `wild_weights` | Ported |
| `Panel.VR.R` | `panel_variance_ratio` | Ported |
| `Spec.shape.R` | `spectral_shape_test` | Ported |
| `Subsample.test.R` | `subsample_variance_ratio` | Ported |
| `VR.minus.1.R` | `variance_ratio_minus_one` | Ported |
| `VR.plot.R`, `stat.plot.R` | `variance_ratio_curve` | Numerical output ported; plotting omitted |
| `WK_stat1.R`, `WK_stat2.R` | internal Whang-Kim routines | Ported |
| `Wald.R`, `Wald1.R`, `covmat.R` | `wald_test` and internal covariance matrix | Ported |
| `Wright.R`, `Wright_stat.R`, `stat.R` | `wright_tests` and internal standardized statistic | Ported |
| `Wright.crit.R` | `wright_critical_values` | Ported |
| `compweexp.R` | internal generalized-spectral weight construction | Ported |
| `kfunc.R` | `quadratic_spectral_kernel` | Ported |

## Exclusions

- R graphics calls are not translated.
- R object naming, row names, and column names are not part of the numerical API.
- The `exrates.rda` data object is retained under `original` but is not converted
  into a compiled Fortran data module.

No exported computational test from `NAMESPACE` is omitted.
