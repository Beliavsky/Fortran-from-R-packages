# betafunctions-fortran

Modern Fortran translation of the computational core of the R package
`betafunctions` 1.9.0 by Haakon Eidem Haakstad.

The upstream package declares `License: CC0`. This port preserves that licensing,
uses SPDX `CC0-1.0` identifiers, and bundles the upstream source tree under
`upstream/betafunctions-master/` for attribution and auditability.

## Scope

Translated numerical functionality includes:

- two- and four-parameter beta density, CDF, quantile, and RNG routines;
- mean/variance parameterization (`AMS`, `BMS`, and related location formulas);
- beta, binomial, beta-binomial, and observed moments;
- two- and four-parameter beta moment fitting;
- beta-binomial density/CDF/RNG;
- generalized real-valued Gamma-Binomial density/CDF/quantile/RNG;
- Lord's two-term compound-binomial approximation and beta-compound-binomial mixtures;
- Livingston-Lewis effective test length and true-score fitting;
- Hanson-Brennan true-score moments and Lord's k;
- binary and multicategory classification accuracy/consistency calculations;
- numerical ROC tables and AUC;
- confusion-matrix statistics;
- Cronbach's alpha and McDonald's omega;
- chi-square model-fit helpers;
- falling/rising factorial and true-score-moment utilities.

Plotting and presentation-only code is intentionally omitted. In particular,
`Beta.gfx.poly.*`, `mdlfit.gfx`, and the output-formatting helper
`MC.out.tabular` are not ported.

## Build

With Fortran Package Manager:

```text
fpm build
fpm test
fpm run --example basic
```

The port has no external package dependencies.

## Main modules

- `betafunctions`: convenience umbrella module.
- `bf_distributions`: probability distributions and RNGs.
- `bf_moments`: moment calculations and beta/true-score fitting.
- `bf_classification`: Livingston-Lewis, Hanson-Brennan, reliability, ROC, and model fit.
- `bf_special`: incomplete beta/gamma and binomial numerical primitives.

Example:

```fortran
program example_beta
  use betafunctions
  implicit none
  real(dp) :: f

  f = beta4_pdf(0.5_dp, 0.0_dp, 1.0_dp, 5.0_dp, 3.0_dp)
  print *, f
end program example_beta
```

## R-to-Fortran naming map

The Fortran API uses descriptive underscore names rather than punctuation-heavy R
names. Important mappings include:

| R | Fortran |
|---|---|
| `dBeta.4P` | `beta4_pdf` |
| `pBeta.4P` | `beta4_cdf` |
| `qBeta.4P` | `beta4_quantile` |
| `rBeta.4P` | `beta4_random` |
| `d/p/q/rBetaMS` | `beta_ms_pdf/cdf/quantile/random` |
| `d/p/rBetaBinom` | `beta_binomial_pmf/cdf/random` |
| `d/p/q/rGammaBinom` | `gamma_binomial_pdf/cdf/quantile/random` |
| `d/p/r cBinom` | `compound_binomial_pmf/cdf/random` |
| `d/rBetacBinom` | `beta_compound_binomial_pmf/random` |
| `Beta.2p.fit` | `beta2_fit` |
| `Beta.4p.fit` | `beta4_fit` |
| `Beta.tp.fit` | `beta_true_score_fit` |
| `HB.beta.tp.fit` | `hb_beta_true_score_fit` |
| `LL.CA`, `LL.CA.MC` | `ll_classify`, `ll_classify_params` |
| `HB.CA`, `HB.CA.MC` | `hb_classify`, `hb_classify_params` |
| `LL.ROC` | `ll_roc`, `ll_roc_params` |
| `HB.ROC` | `hb_roc`, `hb_roc_params` |
| `cba` | `cronbach_alpha` |
| `mdo` | `mcdonald_omega` |
| `Lords.k` | `lords_k` |
| `ETL` | `etl` |
| `R.ETL` | `reliability_from_etl` |
| `dfac`, `afac` | `descending_factorial`, `ascending_factorial` |
| `AUC` | `auc` |

The array-based classification routines handle one or many cut points, so separate
Fortran implementations of the R binary and `*.MC` front ends are unnecessary.

## Validation

The test suite covers independent reference values for special functions and
probability distributions, four-parameter fitting, Livingston-Lewis and
Hanson-Brennan classification matrices, reliability calculations, and RNG behavior.
See `PORTING_NOTES.md` for numerical and compatibility details.
