# survey-fortran

Modern Fortran 2018/FPM translation of the computational core of the R package
`survey` 4.5 (Analysis of Complex Survey Samples).

This port is intentionally array/design-matrix based. It translates numerical
survey methodology rather than R's formula language, S3 methods, plotting,
database backends, or presentation helpers.

## Implemented numerical areas

- multistage stratified/clustered Taylor linearization with finite-population corrections
- lonely-PSU handling and survey design degrees of freedom
- totals, means, ratios, covariance/variance, CDFs, contingency tables, kappa
- replicate-weight variance plus JK1, JKn, BRR/Fay, and bootstrap construction
- calibration/GREG, raking, bounded logit and sinh calibration, post-stratification, weight trimming
- weighted quantiles including upstream `qrule` families and replicate/Taylor uncertainty
- survey GLMs: Gaussian/identity, binomial/logit, Poisson/log with sandwich covariance
- Wald and first-order Rao-Scott contingency-table tests
- PPS Horvitz-Thompson/Yates-Grundy variance support and common joint-inclusion approximations
- weighted 2SLS survey instrumental-variable regression
- generic maximum pseudo-likelihood estimation using `minqa` and `numDeriv`
- generic nonlinear survey regression (`svynls`-style Gauss-Newton/sandwich API)
- weighted Kaplan-Meier, survey Cox regression, log-rank testing, and survey parametric survival regression
- survey t/rank tests, confidence intervals for proportions, contrasts
- Cronbach alpha, weighted correlations, and survey PCA primitives
- special functions used for chi-square/F inference, including weighted-mixture Satterthwaite and saddlepoint tails
- ordinal survey regression (`svyolr` core): logistic, probit, cloglog/Gumbel, and cauchit links
- survey loglinear models with nested deviance/score comparisons
- survey ML factor analysis with effective sample size and varimax
- two-phase/multiphase Dcheck construction, phase variance decomposition, totals/means, and calibration-space projection
- dual-frame/multiframe constant and expected overlap estimators with HT variance
- Preston rescaled multistage bootstrap replicate weights
- pseudo-score, working Rao-Scott score, Wald term, and misspecification-adjusted LRT model tests

See `docs/TRANSLATION_COVERAGE.md` for detailed coverage and known gaps.

## Build

With FPM installed:

```sh
fpm test
fpm run --example basic_survey
```

The project has also been validated directly with GNU Fortran using strict
Fortran 2018 checks; see `scripts/strict_test.sh`.

## Minimal example

```fortran
program basic_survey
    use survey
    implicit none

    type(survey_design_t) :: d
    type(svystat_t) :: m
    real(dp) :: y(4,1), w(4)
    integer :: psu(4,1), strata(4,1)

    y(:,1) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    w = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp]
    psu(:,1) = [1, 2, 3, 4]
    strata(:,1) = 1

    call make_design(w, psu, d, strata=strata)
    m = svy_mean(y, d)

    print '(a,f12.6)', 'weighted mean = ', m%estimate(1)
    print '(a,f12.6)', 'SE            = ', sqrt(m%variance(1,1))
end program basic_survey
```

## Design philosophy

R objects such as formulas, model frames, factors, `Surv` objects, S3/S4
classes, and DBI-backed designs are represented by explicit numeric arrays and
typed Fortran results. This keeps the statistical core usable in standalone
Fortran applications and makes ownership of allocation and model matrices
explicit.

## Licensing

The upstream `survey` package is GPL-2 or GPL-3. This v0.2.0 release preserves that licensing choice. Original attribution is in
`LICENSE`; full texts and dependency provenance are in `LICENSES.md` and
`licenses/`.
