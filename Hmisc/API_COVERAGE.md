# API coverage

This translation concentrates on portable computational code and intentionally omits plotting, printing, formatting, HTML/LaTeX/Typst generation, data download/import interfaces, interactive display, and R-specific object-system glue.

The current package-wide coverage basis contains **142 exported computational R functions**. **136 are mapped (95.8%)** and **6 remain untranslated**. The audited denominator excludes presentation-only `ordGridFun` and `nFm`, plus `var.inner`, which is an internal formula-parsing helper rather than a standalone computational API. The exact untranslated names and one mapping record for every translated R function are recorded in `fpm.toml`; README.md carries the matching summary table. The fraction measures mapped computational functions rather than full R compatibility.

Mapped numerical areas now include:

- weighted means, variances, frequency tables, quantiles, and ranks;
- Pearson/Spearman correlations and Hoeffding dependence statistics;
- censored and paired censored concordance;
- Somers' D;
- Gini mean difference and pseudomedian;
- grouping and nearest-value helpers;
- trapezoidal integration and binary sample-size calculation;
- dual standard deviations;
- binary power, sample-size, and allocation formulas;
- two-point Weibull, Gompertz, and lognormal survival fitting;
- the event, non-event, and total NRI core of `improveProb`;
- caller-supplied-knot restricted cubic spline basis evaluation;
- proportional-odds power, sample size, and category-probability transformation;
- grouped means, standalone Spearman correlation, step-function evaluation, and numeric sorting/deduplication;
- coincident-point counting and matrix-vector prediction helpers.
- survival-study power calculations for one treatment contrast and two-stratum interactions;
- weighted empirical CDFs and mean/SD/median interval summaries;
- linear interpolation with endpoint extrapolation;
- James-Stein shrinkage of grouped means;
- Fleiss-Tytun-Ury binary power and sample-size approximations;
- empirical-CDF coordinate construction, numeric lagging, and missing-value filtering/filling;
- cumulative-category indicator encoding.
- cluster design effects and intracluster correlation estimation;
- generalized likelihood-based R-squared measures with effective sample sizes;
- Wilson/asymptotic binomial confidence intervals, linear tabulated inversion, and binary-item scoring cores;
- Mantel-Haenszel risk ratios and cumulative diagnostic likelihood ratios;
- numerator/denominator expansion, multinomial sampling with explicit uniform draws, and centered QR transforms;
- QR least-squares/inversion, first-principal-component scoring, normal-theory confidence summaries, and deterministic k-nearest selection;
- bare linear-model fitting, absolute prediction-error diagnostics, and numeric censored-concordance summaries;
- two-group log-rank statistics, numeric Spearman rho-squared tests, display-bin width redistribution, and deterministic numeric tie jittering;
- Harrell-Davis quantiles, Pearson contingency-table chi-square testing, numeric `conTestkw`, and clustered two-sample mean testing.
- right-censored Kaplan-Meier curves and step evaluation, default one- and two-way grouped means, and character digit/numeric tests.
- tolerance-based numeric matching, hierarchical sequential-condition assignment, and paired-difference confidence calculations;
- bootstrap percentile confidence limits for means and tabulated-inverse smearing estimates;
- Bezier-curve evaluation, explicit-cut categorization, deterministic proportional-odds cut simulation, and restricted-cubic-spline function/restate numerical cores.
- numeric `spearman.test`, integer-coded category pooling and Pearson chi-square summaries, and regular complete-case principal-components analysis.
- deterministic Kaplan-Meier bootstrap summaries, scalar case-control matching, and rexhaustive largest-empty-rectangle geometry.
- Gaussian conjugate Bayesian updating, Gaussian-mixture predictive/posterior density and CDF evaluation, posterior mixture means, and analytic single-prior power, and deterministic two-component-mixture posterior-tail power.

- tabulated-prior `gbayes2` Bayesian integration, proportional-odds category projections from counts, adjacent-time transition count/probability kernels, and numeric 0/1 `ynbind` binding/sorting.
- missingness-driven variable ordering and imputation-count heuristics, grouped matrix means, and exact proportional-odds Markov occupancy propagation with absorbing states.
- robust LOWESS smoothing, weighted local-polynomial smoothing, ordinal bootstrap grouping with approximate or supplied-resample coverage checks, observed-point smoothing of grouped curves, and deterministic ordinal Markov simulation with supplied uniforms.
- numeric median/constant imputation, default raw fixed-sample moving means/quartiles, and multi-key grouped arithmetic means.
- linear/spline/categorical `aregTran` design transforms, the all-linear numeric `areg` fit path, numeric data-frame reduction and representativeness counts, variable-clustering similarities/linkage, linear redundancy elimination, sequential Gaussian posterior assertion probabilities, and model-agnostic Markov transition occupancy propagation.
- linear `areg.boot` optimism correction, Rubin multiple-imputation combination, numeric stored-imputation insertion, sequential two-group mean-difference estimation, and the linear-time fixed-x residual-bootstrap core of `rm.boot`.
- numeric proportional-odds likelihood-ratio testing, intercept-only Markov occupancy calibration, and supplied-replicate ordinal-regression power calculations.

Material compatibility gaps inside mapped functions are documented in each `[[extra.translation.function]]` record. Examples include omitted R object construction and P-value matrices, alternate weighted-quantile interpolation types, large-sample randomized/bootstrap behavior in `pMedian`, and the additional uncertainty/IDI outputs of `improveProb`.

The package is a **substantial** computational translation, not a complete R-interface reimplementation. The remaining untranslated exports are now concentrated in the heaviest model-fitting, imputation, clustering, sequential-Markov, and transformation frameworks: `aregImpute`, `biVar`, `curveRep`, `estSeqMarkovOrd`, `transace`, and `transcan`.