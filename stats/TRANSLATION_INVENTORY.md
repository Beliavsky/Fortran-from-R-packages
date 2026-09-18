# R `stats` translation inventory

This first-pass inventory records numerical functionality already present in
the MIT-licensed `R-to-Fortran` runtime at revision
`6c930904ea2517703550ff73ebb5a4f3c6b867aa`. It guides extraction; presence in
the source runtime does not by itself establish full parity with R.

| Area | Representative existing procedures | Intended disposition |
|---|---|---|
| Sample quantile estimators | `quantile`, `IQR` | All nine unweighted Hyndman-Fan types now delegate to pure shared `r_quantile`, with deterministic R fixtures; type 7 remains the default. |
| Serial-correlation diagnostics | `Box.test` | Pure Box-Pierce and Ljung-Box statistics reuse the existing ACF and support fitted-parameter adjustments, with deterministic R fixtures at several lags. |
| One-way mean inference | `oneway.test` | Welch unequal-variance tests and pooled-variance ANOVA use explicit integer groups and retain fractional degrees of freedom. Deterministic R fixtures cover two and three groups. |
| Robust scale testing | `fligner.test` | Median-centered normal rank scores reuse the shared median, ranking, normal quantile, and chi-square tail functions. R fixtures cover ties, unequal group sizes, and singleton groups. |
| Paired categorical inference | `mcnemar.test` | Pure matrix interfaces support McNemar correction and the larger-square-table symmetry test, with deterministic R fixtures. |
| Multiple testing | `p.adjust` | Pure Holm, Hochberg, Hommel, Bonferroni, BH/fdr, and BY corrections reuse shared ordering; ties, missing values, and larger hypothesis counts have R fixtures. All method names are supported. |
| Descriptive statistics | `median`, `quantile`, `sd`, `var`, `cov`, `cor`, `cov2cor`, `IQR`, `mad`, `weighted.mean` | Pure user-facing scalar, vector, and matrix adapters reuse the canonical `rfortran-core` implementations and have deterministic R fixtures. |
| Empirical distributions | `ecdf`, `density`, `bw.nrd0`, `bw.nrd` | Compact ECDF objects, normal-reference bandwidths, and direct weighted kernel densities have deterministic R comparisons. A plotting-free histogram companion is provided but is not counted because R's `hist` belongs to `graphics`. |
| Interpolation and monotone fitting | `approx`, `approxfun`, `isoreg`, `spline`, `splinefun` | Pure linear, constant, FMM cubic, and natural cubic interpolation objects handle sorting, averaged ties, endpoint rules, derivatives, and extrapolation; weighted PAVA supplies monotone least-squares fits. |
| Probability distributions | Normal, t, F, chi-square, beta, gamma, count and rank distributions | Reuse `rfortran-core`; do not copy runtime implementations. |
| Linear models and ANOVA | `lm_fit_general`, prediction, confidence intervals, Cook's distance, AIC, one-way ANOVA, nested-model comparison | Dense weighted and rank-deficient fitting is implemented with alias metadata. Formula-free one-way and nested-model F tests have deterministic R fixtures. |
| Nonlinear least squares | Formula-free model callbacks, Gauss-Newton iteration, bounds, covariance | `nls_fit` implements analytic or finite-difference Jacobians, active bounds, rank diagnostics, and deterministic R fixtures; `nls_fit_plinear` profiles linear coefficients, `nls_confint_profile` adds profile-RSS intervals, and five common self-start curve families have pure value and gradient helpers. |
| Generalized linear models | Six GLM families, eight links, prediction and residuals | A shared dense IRLS engine handles binomial, Poisson, Gaussian, Gamma and their two quasi variants, prior weights, offsets, optional intercepts, rank-deficient designs, dispersion, deviance, covariance, and deterministic R fixtures. |
| Classical tests | t, chi-square, proportion, correlation, exact binomial and Poisson, variance ratio, Bartlett, Friedman, Fisher, Wilcoxon, Kruskal-Wallis and KS tests | All listed families have typed initial implementations and deterministic R fixtures; rank-based grouped procedures include tie correction. |
| Multivariate methods | PCA, Mahalanobis distance, scaling, multivariate regression, MANOVA, canonical correlation | Dense full-rank multivariate response fitting, all four standard MANOVA criteria, and canonical correlations now accompany the existing PCA, distance, and scaling procedures. |
| Clustering | `kmeans`, `hclust`, `cutree` and distance calculations | Implemented with deterministic Lloyd fitting and adapters around canonical `fastcluster` distances and linkage. |
| Smoothing | running medians, kernel smoothing, LOWESS, LOESS and smoothing splines | `runmed`, `ksmooth`, `lowess`, univariate direct-surface LOESS, and natural cubic smoothing splines are implemented with R fixtures. Smoothing splines include reduced-knot bases, duplicate aggregation, and automatic GCV selection. |
| Time series | ACF/CCF, AR fitting, ARIMA simulation/fitting/prediction, filters and decomposition | Filtering, univariate ACF/CCF/PACF, theoretical ARMA ACF, deterministic nonseasonal integrated simulation, AR fitting and spectra, periodograms, classical and STL decomposition, Holt-Winters, and multiplicative seasonal ARIMA with CSS or exact Gaussian likelihood, stationary and integrated response gaps, regressors, fixed coefficients, coefficient covariance, point forecasts, standard errors, and normal intervals are implemented. |
| Optimization | `uniroot`, `optimize`; multidimensional BFGS, L-BFGS-B, conjugate gradient, Nelder-Mead, simulated annealing and constrained optimization | Pure bracketed root finding with optional interval expansion and scalar safeguarded Brent optimization are implemented with typed diagnostics and deterministic R fixtures. Multidimensional methods should use canonical optimization packages through future `stats`-style adapters. |
| Integration | `integrate` | Pure adaptive Simpson quadrature supports finite, semi-infinite, and whole-line bounds, explicit error tolerances, subdivision limits, and typed diagnostics. Vectorized integrands and QUADPACK-specific behavior remain omitted. |
| R infrastructure | formulas, model frames, terms, S3 dispatch, printing and plotting | Omit initially or replace with explicit typed Fortran APIs. |

## Existing validation assets

`R-to-Fortran/test` contains focused programs for linear models, generalized
linear models, diagnostics, inferential statistics, clustering, PCA,
optimization, integration, smoothing, and time-series models. Tests should be
adapted with each extracted slice, and independent comparisons against R
should be added before a function is described as complete.

## Next extraction

Seasonal simulation would be a Fortran extension beyond R's nonseasonal
`arima.sim`. Formula processing remains a cross-cutting omission. The next GLM
work is binomial count-matrix input, additional families, and user-controlled
starts and convergence settings. Separable nonlinear models
and five common self-start curve families are implemented; automatic
starting-value estimation remains out of scope.
It should reuse canonical translated engines where licensing and dependency
direction permit, rather than copying their implementations.
