# Translation coverage

Status: **partial initial implementation**.

The manifest uses 465 public namespace exports from R 4.6.0 as its denominator.
Partial mappings count as represented, so the fraction measures API coverage,
not the fraction of all behavior implemented. Some method-level counterparts
are documented below without being counted as separate primary functions.

| R function | Fortran counterpart | Status | Notes |
|---|---|---|---|
| `Box.test` | `box_test` | Substantial | Box-Pierce (default) and Ljung-Box statistics support fitted-parameter degrees-of-freedom adjustments. Requires finite input and `0 <= fitdf < lag < size(x)`. Missing-value handling and partial method names are omitted. |
| `oneway.test` | `oneway_test` | Substantial | Welch unequal-variance ANOVA (default) or pooled-variance inference accepts finite observations and arbitrary integer labels, with at least two observations per group. Real degrees of freedom and input/variance status are retained. Formula processing and missing-value omission are omitted. |
| `fligner.test` | `fligner_test` | Substantial | Median-centered Fligner-Killeen testing uses normal scores of average ranks. Finite samples support arbitrary integer group labels, unequal sizes, and singleton groups. Formula/list input and missing-value omission are omitted. |
| `mcnemar.test` | `mcnemar_test` | Substantial | Integer or real square tables support the two-category continuity correction and larger-table symmetry statistic. Raw paired observations and factor conversion are omitted. Zero opposite-cell totals return NaNs, as in R. |
| `p.adjust` | `p_adjust` | Substantial | All methods (Holm, Hochberg, Hommel, Bonferroni, BH/fdr, BY, and none) preserve ties, input order, NaNs, and explicit hypothesis counts. Names and partial method-name matching are omitted; finite inputs must be in `[0,1]`. |
| `median` | `median` | Substantial | Pure real-vector median with optional NaN removal. |
| `quantile` | `quantile` | Partial | Scalar or vector probabilities support all nine Hyndman-Fan estimators, defaulting to type 7. Names, factors, configurable fuzz, and R's near-boundary probability tolerance are omitted. Empty samples and unremoved NaNs return NaN. |
| `sd` | `sd` | Substantial | Pure real-vector sample standard deviation with optional NaN removal. |
| `var` | `var` | Substantial | Vector sample variance or covariance matrix of data columns. |
| `cov` | `cov` | Partial | Vector covariance or column covariance matrix; optional removal uses pairwise complete observations. Rank methods are omitted. |
| `cor` | `cor` | Partial | Vector Pearson correlation or column correlation matrix; Kendall and Spearman methods are omitted here. |
| `cov2cor` | `cov2cor` | Substantial | Pure conversion of a square covariance matrix. |
| `IQR` | `iqr` | Substantial | Interquartile range supports all nine quantile estimators, defaulting to type 7, with optional NaN removal. |
| `mad` | `mad` | Partial | Configurable center and scale constant; the low/high even-sample variants are omitted. |
| `weighted.mean` | `weighted_mean` | Substantial | Pure real-vector weighted mean with optional NaN removal. |
| `ecdf` | `ecdf`, `predict_ecdf` | Substantial | Compact sorted ECDF objects preserve duplicate frequencies and support vector evaluation. R step-function closures and plotting methods are omitted. |
| `bw.nrd0` | `bw_nrd0` | Substantial | R's robust default normal-reference rule and zero-spread fallbacks are translated. |
| `bw.nrd` | `bw_nrd` | Substantial | R's normal-reference bandwidth rule is translated. |
| `density` | `density` | Partial | Direct weighted KDE on a regular grid supports Gaussian, rectangular, triangular, and Epanechnikov kernels. R's FFT implementation, automatic missing-value handling, and additional kernels and bandwidth selectors are omitted. |
| `approx` | `approx` | Substantial | Explicit output points, linear or constant interpolation, mean resolution of duplicate predictors, separate endpoint rules, and the constant mixing fraction are supported. Missing-value omission and generated default grids are omitted. |
| `approxfun` | `approxfun`, `predict_approxfun` | Substantial | A reusable derived type replaces R's closure and supports the translated `approx` controls. |
| `isoreg` | `isoreg` | Substantial | Pure increasing PAVA fits return fitted levels and block endpoints for strictly increasing predictors; positive observation weights are a Fortran extension. |
| `spline` | `spline` | Partial | Explicit output points use FMM or natural cubic interpolation with averaged ties and derivatives through order three. Generated default grids, periodic splines, and Hyman filtering are omitted. |
| `splinefun` | `splinefun`, `predict_splinefun` | Partial | A reusable derived type replaces R's closure and reproduces FMM cubic or natural linear extrapolation. Periodic and monotone Hermite methods are omitted. |
| `uniroot` | `uniroot` | Partial | Pure formula-free bracketed bisection returns the root, function value, precision, and convergence diagnostics. Optional symmetric, lower, or upper interval expansion is supported; directional monotonicity hints, supplied endpoint values, tracing, and callback data are omitted. |
| `optimize` | `optimize` | Substantial | Safeguarded Brent search minimizes or maximizes a pure scalar callback over a finite interval and returns typed diagnostics. Additional callback arguments are omitted. |
| `integrate` | `integrate` | Substantial | Adaptive Simpson quadrature handles pure scalar callbacks, finite bounds, either semi-infinite interval, or the whole real line with absolute and relative tolerances. Vectorized callback evaluation and QUADPACK-specific diagnostics are omitted. |
| `lm.fit`, `lm.wfit` | `lm_fit`, `lm_fit_general` | Substantial | Dense real least squares with an optional intercept, nonnegative prior weights, numerical rank detection, deterministic alias selection, and minimum-subspace fitting. |
| `aov` | `one_way_anova` | Partial | Formula-free one-way ANOVA for explicit contiguous integer groups. Contrasts, multiple strata, and multivariate responses are omitted. |
| `anova.lm` | `anova_lm` | Partial | Extra-sum-of-squares F comparison for two nested fits using identical responses and weights. Sequential term tables and formula dispatch are omitted. |
| `predict.lm` | `predict_lm`, `lm_predict_general` | Partial | Point predictions only. |
| `nls` | `nls_fit`, `nls_fit_plinear` | Partial | Pure formula-free Gauss-Newton fitting accepts model callbacks, optional analytic Jacobians, bounds, active-set projection, step halving, finite-difference fallback, rank diagnostics, covariance, and standard errors. Separable models profile their linear coefficients. Formula processing and the complete `port` algorithm are omitted. |
| `confint.nls` | `nls_confint_profile` | Partial | Pure profile-RSS intervals refit the free parameters and use a Student-t threshold. Formula methods and R profile objects are omitted. |
| `SSasymp` | `ss_asymp`, `ss_asymp_model`, `ss_asymp_jacobian` | Substantial | Pure elemental curve and NLS callbacks; automatic starting-value estimation is omitted. |
| `SSgompertz` | `ss_gompertz`, `ss_gompertz_model`, `ss_gompertz_jacobian` | Substantial | Pure elemental curve and NLS callbacks; automatic starting-value estimation is omitted. |
| `SSlogis` | `ss_logis`, `ss_logis_model`, `ss_logis_jacobian` | Substantial | Pure elemental curve and NLS callbacks; automatic starting-value estimation is omitted. |
| `SSmicmen` | `ss_micmen`, `ss_micmen_model`, `ss_micmen_jacobian` | Substantial | Pure elemental curve and NLS callbacks; automatic starting-value estimation is omitted. |
| `SSweibull` | `ss_weibull`, `ss_weibull_model`, `ss_weibull_jacobian` | Substantial | Pure elemental curve and NLS callbacks; automatic starting-value estimation is omitted. |
| `coef.lm` | `lm_coef`, `lm_fit_t%coef` | Substantial | Numeric coefficient access without names or alias metadata. |
| `summary.lm` | fields of `lm_fit_t` | Partial | Fitted values, residuals, sigma, degrees of freedom, R-squared, and adjusted R-squared. |
| `confint.lm` | `lm_confint` | Partial | Coefficient intervals for estimable coefficients; aliased coefficients return NaN intervals. |
| `cooks.distance.lm` | `lm_cooks_distance` | Partial | Ordinary and prior-weighted dense linear models. |
| `extractAIC.lm` | `lm_aic` | Partial | Constant-free Gaussian score suitable for comparing fits on identical data. |
| `glm.fit` | `glm_fit`, family-specific fitters | Substantial | Dense IRLS for binomial, Poisson, Gaussian, Gamma, quasi-binomial, and quasi-Poisson families; prior weights, offsets, optional intercepts, rank detection, aliases, covariance, dispersion, and deviance are retained. Starts, control objects, missing values, and arbitrary family callbacks are omitted. |
| `predict.glm` | `glm_predict_response` | Partial | Response-scale predictions for all supported families and links. |
| `residuals.glm` | `glm_pearson_resid` | Partial | Prior-weighted Pearson residuals for all supported families. |
| `binomial` | embedded in `glm_binomial_fit` | Substantial | Logit, probit, complementary-log-log, and cauchit links with binomial variance, but no dynamic family object or two-column count response. |
| `poisson` | embedded in `glm_poisson_fit` | Substantial | Log, identity, and square-root links with Poisson variance, but no dynamic family object. |
| `gaussian` | embedded in `glm_gaussian_fit` | Partial | Identity, log, and inverse links with estimated dispersion, without a dynamic family object. |
| `Gamma` | embedded in `glm_gamma_fit` | Partial | Inverse, log, and identity links with estimated dispersion, without a dynamic family object. |
| `quasibinomial` | `glm_quasibinomial_fit` | Substantial | Supported binomial links, Pearson dispersion, and scaled covariance without a dynamic family object. |
| `quasipoisson` | `glm_quasipoisson_fit` | Substantial | Supported Poisson links, Pearson dispersion, and scaled covariance without a dynamic family object. |
| `make.link` | `glm_link_value`, `glm_inverse_link`, `glm_mu_eta` | Partial | Pure elemental helpers implement identity, log, inverse, logit, probit, complementary-log-log, cauchit, and square-root links. |
| `t.test` | `t_test`, `t_test_one`, `t_test_two`, `t_test_p_value` | Partial | Two-sided one-sample, Welch, pooled-variance, and paired tests with 95 percent intervals. |
| `chisq.test` | `chisq_test` | Partial | Goodness-of-fit and contingency-table tests, including Yates correction. |
| `prop.test` | `prop_test` | Partial | One- and multiple-sample two-sided tests with score confidence intervals. |
| `cor.test` | `cor_test` | Partial | Pearson tests with Fisher intervals and asymptotic Spearman tests. |
| `binom.test` | `binom_test` | Partial | Exact two-sided binomial probabilities and Clopper-Pearson confidence intervals are supported. Directional alternatives are omitted. |
| `poisson.test` | `poisson_test` | Partial | Exact two-sided one-sample rate inference and chi-square confidence intervals are supported. Multiple-sample conditional tests and directional alternatives are omitted. |
| `var.test` | `var_test` | Partial | Two-sided F inference for two explicit finite samples supports configurable null ratios and confidence levels. Formula input and directional alternatives are omitted. |
| `bartlett.test` | `bartlett_test` | Substantial | Bartlett's chi-square test supports finite observations and explicit positive integer group labels. Formula and list inputs are omitted. |
| `friedman.test` | `friedman_test` | Substantial | Complete block-by-treatment matrices are tested with average ranks and tie correction. Formula and unstacked-vector inputs are omitted. |
| `fisher.test` | `fisher_test` | Partial | Two-sided exact inference and conditional odds estimate for 2-by-2 tables. |
| `wilcox.test` | `wilcox_test` | Partial | Asymptotic rank-sum and paired signed-rank tests with tie corrections. |
| `kruskal.test` | `kruskal_test` | Substantial | Tie-corrected test for explicit values and integer group labels. |
| `ks.test` | `ks_test`, `ks_test_normal`, `ks_test_two_sample` | Partial | Asymptotic normal-reference and two-sample tests. |
| `prcomp` | `prcomp` | Substantial | Dense real SVD-based PCA with centering, scaling, rank truncation, and scores. |
| `scale` | `scale`, `scale_columns` | Partial | Logical column centering and scaling, with transformation metadata. |
| `mahalanobis` | `mahalanobis` | Substantial | Squared distances from covariance or inverse-covariance matrices. |
| multivariate `lm` | `mlm_fit` | Partial | Full-rank dense response matrices with an optional intercept, fitted values, residuals, and residual SSCP. Weights and rank-deficient designs are omitted. |
| `manova`, `summary.manova` | `manova_lm` | Partial | Two nested multivariate fits produce Pillai, Wilks, Hotelling-Lawley, and Roy criteria, R-compatible approximate F statistics, degrees of freedom, and p-values. Formula terms and sequential tables are omitted. |
| `cancor` | `cancor` | Substantial | Full-rank dense matrices, optional centering, canonical correlations, coefficient matrices, and retained centers. Explicit center vectors and rank-deficient inputs are omitted. |
| `dist` | `dist`, `dist_matrix` | Substantial | Six metrics with R-layout condensed or full-matrix results. |
| `kmeans` | `kmeans` | Partial | Deterministic Lloyd fitting from explicit initial centers. |
| `hclust` | `hclust`, `hclust_condensed` | Substantial | Shared `fastcluster` implementation with R merge conventions. |
| `cutree` | `cutree` | Partial | Scalar cluster-count or height cuts. |
| `runmed` | `runmed` | Substantial | Odd windows and the `median`, `keep`, and `constant` endpoint rules. |
| `ksmooth` | `ksmooth` | Substantial | Normal and box kernels at explicit or default evaluation locations. |
| `lowess` | `lowess` | Substantial | Robust locally linear smoothing with configurable span, iterations, and delta. |
| `loess`, `predict.loess` | `loess_fit`, `predict_loess` | Substantial | Pure univariate direct-surface local polynomials of degree zero through two with tricube neighborhoods, prior weights, Gaussian or Tukey-bisquare robust fitting, retained diagnostics, and prediction at arbitrary finite locations. Formula input, multiple predictors, interpolation surfaces, standard errors, and influence diagnostics are omitted. |
| `smooth.spline`, `predict.smooth.spline` | `smooth_spline`, `predict_smooth_spline` | Substantial | Pure natural cubic smoothing splines support R-compatible default reduced-knot counts, explicit reduced or all-knots bases, prior weights, tolerance-based duplicate aggregation, automatic GCV or ordinary-CV selection, and explicit `spar`, positive `lambda`, or target degrees of freedom. Predictions cover derivative orders zero through two with natural extrapolation. Formula methods, arbitrary numeric knot vectors, and standard-error calculations are omitted. |
| `filter` | `filter`, `filter_linear`, `filter_recursive` | Substantial | One-sided and centered convolution plus recursive filtering with optional initial states; circular filtering is omitted. |
| `acf` | `acf` | Partial | Univariate correlations or biased covariances with explicit lag metadata. |
| `pacf` | `pacf` | Partial | Univariate Yule-Walker MLE estimates through the requested lag. |
| `ccf` | `ccf` | Partial | Univariate cross-correlations or biased cross-covariances at positive and negative lags. |
| `ARMAacf` | `arma_acf` | Substantial | Theoretical stationary ARMA autocorrelations from finite covariance equations. |
| `arima.sim` | `arima_sim` | Substantial | Nonseasonal ARIMA simulation from explicit retained and burn-in innovations, including R-compatible ordinary integration. RNG callbacks, automatic burn-in selection, and time-series attributes are omitted. |
| `arima` | `arima_css`, `arima_ml` | Substantial | Multiplicative seasonal CSS and exact Gaussian innovations fits with ordinary and seasonal integration, an optional stationary mean, user regressors, explicit fixed-parameter values and masks, and observed numerical-Hessian covariance. Exact ML handles internal response gaps, using Gaussian conditioning for stationary models and exact diffuse state-space initialization for integrated models. |
| `predict.Arima` | `predict_arima`, `predict_arima_interval` | Substantial | Finite-history exact and conditional CSS point forecasts accept future regressors and ordinary and seasonal integration; innovation-based standard errors and configurable normal intervals are also returned. Coefficient-estimation uncertainty is omitted. |
| `ar` | `ar`, `ar_yw`, `ar_burg`, `ar_ols`, `ar_mle` | Substantial | Univariate Yule-Walker, Burg, unconstrained OLS, and exact stationary Gaussian ML fitting with fixed-order or AIC selection and typed fit metadata; multivariate fitting is omitted, and ML defaults to exact rather than `ar.mle`'s historical approximate likelihood. |
| `spec.ar` | `spec_ar` | Substantial | Univariate spectra from fixed-order or AIC-selected Yule-Walker, Burg, OLS, or ML fits; multivariate fits are omitted. |
| `spec.pgram` | `spec_pgram` | Substantial | Raw and repeated modified-Daniell-smoothed univariate or multivariate periodograms, including cross-spectra, coherence, phase, tapering, padding, scaling, and metadata. |
| `spectrum` | `spectrum` | Partial | Default-method adapter for the translated `spec_pgram` path. |
| `decompose` | `decompose` | Substantial | Classical additive or multiplicative decomposition using R's default moving average or caller-supplied centered convolution weights. |
| `stl` | `stl` | Substantial | Seasonal, trend, and low-pass LOESS controls, periodic seasonality, jump interpolation, inner iterations, and robust bisquare reweighting. |
| `HoltWinters` | `holt_winters`, `holt_winters_auto`, `holt_winters_optimize`, `holt_winters_reduced` | Substantial | Additive and multiplicative seasonal recurrences, automatic initialization, bounded estimation, Holt linear trend, simple exponential smoothing, and seasonal models without trend. |
| `predict.HoltWinters` | `predict_holt_winters`, `predict_holt_winters_interval` | Substantial | Point forecasts and configurable normal prediction intervals using R's variance formulas; multiplicative intervals require trend. |

The manifest counts eighty-nine primary R `stats` functions represented by the current descriptive-statistics,
linear-model, generalized-linear-model, testing, multivariate, clustering, smoothing, and time-series slices. The supporting `scale`
translation belongs to R's `base` package and is not included in that count. Linear-model accessors and diagnostics
will be counted separately after their behavioral parity tests cover the
corresponding R APIs.
