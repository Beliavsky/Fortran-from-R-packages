! SPDX-License-Identifier: MIT
! SPDX-FileComment: Result types for the modern Fortran translation of R stats functionality.
module r_stats_types
   use r_kinds, only: dp
   implicit none
   private

   public :: anova_comparison_t, chisq_test_result_t, cor_test_result_t
   public :: density_result_t, ecdf_t, fisher_test_result_t, glm_fit_t, histogram_t
   public :: exact_count_test_result_t, variance_test_result_t
   public :: cubic_spline_t, interpolation_t, isoreg_t
   public :: holt_winters_fit_t, holt_winters_forecast_t, kmeans_result_t
   public :: kruskal_test_result_t
   public :: ks_test_result_t, lm_fit_t, nls_fit_t, one_way_anova_t, prcomp_fit_t
   public :: acf_result_t, ar_fit_t, arima_fit_t, arima_forecast_t
   public :: multivariate_spectrum_result_t, prop_test_result_t
   public :: cancor_result_t, manova_result_t, mlm_fit_t
   public :: loess_fit_t, scale_result_t, smooth_spline_fit_t, smooth_xy_t
   public :: integrate_result_t, optimize_result_t, uniroot_result_t
   public :: spectrum_result_t
   public :: seasonal_decomposition_t
   public :: stl_result_t
   public :: t_test_result_t
   public :: wilcox_test_result_t

   type :: exact_count_test_result_t
      !! Stores an exact one-sample binomial or Poisson test.
      real(dp) :: statistic = 0.0_dp !! Observed success or event count.
      real(dp) :: parameter = 0.0_dp !! Trial count or exposure.
      real(dp) :: p_value = 1.0_dp !! Exact two-sided probability.
      real(dp) :: estimate = 0.0_dp !! Estimated probability or event rate.
      real(dp) :: null_value = 0.0_dp !! Null probability or event rate.
      real(dp) :: conf_low = 0.0_dp !! Lower exact confidence limit.
      real(dp) :: conf_high = 0.0_dp !! Upper exact confidence limit.
      integer :: method = 0 !! One for binomial and two for Poisson.
      integer :: status = 0 !! Zero for valid input and one otherwise.
   end type exact_count_test_result_t

   type :: variance_test_result_t
      !! Stores a two-sample F test for equality of variances.
      real(dp) :: statistic = 0.0_dp !! Ratio of the two sample variances.
      real(dp) :: degrees_freedom1 = 0.0_dp !! Numerator degrees of freedom.
      real(dp) :: degrees_freedom2 = 0.0_dp !! Denominator degrees of freedom.
      real(dp) :: p_value = 1.0_dp !! Two-sided F-test probability.
      real(dp) :: estimate = 0.0_dp !! Estimated variance ratio.
      real(dp) :: null_value = 1.0_dp !! Null variance ratio.
      real(dp) :: conf_low = 0.0_dp !! Lower confidence limit for the variance ratio.
      real(dp) :: conf_high = 0.0_dp !! Upper confidence limit for the variance ratio.
      integer :: status = 0 !! Zero for valid input and one otherwise.
   end type variance_test_result_t

   type :: uniroot_result_t
      !! Stores a bracketed scalar root search corresponding to R `uniroot`.
      real(dp) :: root = 0.0_dp !! Estimated zero of the objective function.
      real(dp) :: function_value = 0.0_dp !! Objective value at `root`.
      real(dp) :: estimated_precision = 0.0_dp !! Width of the final bracketing interval.
      integer :: iterations = 0 !! Number of bisection refinements performed.
      integer :: status = 0 !! Zero on convergence or a numerical-method status code.
      logical :: converged = .false. !! Whether the requested tolerance was attained.
   end type uniroot_result_t

   type :: optimize_result_t
      !! Stores a bounded scalar optimization result corresponding to R `optimize`.
      real(dp) :: location = 0.0_dp !! Estimated minimizing or maximizing argument.
      real(dp) :: objective = 0.0_dp !! Objective value at `location`.
      integer :: evaluations = 0 !! Number of objective evaluations.
      integer :: status = 0 !! Zero on convergence or a numerical-method status code.
      logical :: converged = .false. !! Whether the requested tolerance was attained.
   end type optimize_result_t

   type :: integrate_result_t
      !! Stores an adaptive finite-interval quadrature result corresponding to R `integrate`.
      real(dp) :: value = 0.0_dp !! Estimated integral.
      real(dp) :: absolute_error = 0.0_dp !! Estimated absolute quadrature error.
      integer :: subdivisions = 0 !! Number of accepted leaf intervals.
      integer :: evaluations = 0 !! Number of integrand evaluations.
      integer :: status = 0 !! Zero on convergence or a numerical-method status code.
      logical :: converged = .false. !! Whether every leaf met its local error target.
   end type integrate_result_t

   type :: acf_result_t
      !! Stores a univariate autocorrelation, autocovariance, or cross-correlation sequence.
      real(dp), allocatable :: lag(:) !! Integer-valued lags represented in time units.
      real(dp), allocatable :: value(:) !! Statistic corresponding to each lag.
      integer :: n_used = 0 !! Number of paired observations used.
      integer :: status = 0 !! Shared-core status; zero indicates success.
      logical :: covariance = .false. !! True when `value` contains covariances rather than correlations.
   end type acf_result_t

   type :: ar_fit_t
      !! Stores a univariate Yule-Walker autoregressive fit.
      real(dp), allocatable :: coefficients(:) !! AR coefficients in increasing lag order.
      real(dp), allocatable :: residuals(:) !! Residuals, with NaNs before the fitted order.
      real(dp), allocatable :: aic(:) !! Relative AIC values for `aic_order`.
      integer, allocatable :: aic_order(:) !! Candidate orders corresponding to `aic`.
      real(dp) :: mean = 0.0_dp !! Fitted series mean.
      real(dp) :: intercept = 0.0_dp !! OLS intercept after optional demeaning.
      real(dp) :: variance = 0.0_dp !! Bias-adjusted innovation variance.
      real(dp) :: log_likelihood = 0.0_dp !! Maximized Gaussian log likelihood when available.
      integer :: order = 0 !! Selected autoregressive order.
      integer :: n_used = 0 !! Number of fitted observations.
      integer :: status = 0 !! Zero on success or a shared-core failure status.
      logical :: demeaned = .false. !! Whether the sample mean was removed.
      character(len=16) :: method = "yule-walker" !! Autoregressive fitting method.
   end type ar_fit_t

   type :: arima_fit_t
      !! Stores a multiplicative seasonal CSS or exact-likelihood ARIMA fit.
      real(dp), allocatable :: ar(:) !! Autoregressive coefficients in increasing lag order.
      real(dp), allocatable :: ma(:) !! Moving-average coefficients in increasing lag order.
      real(dp), allocatable :: seasonal_ar(:) !! Seasonal AR coefficients by seasonal lag.
      real(dp), allocatable :: seasonal_ma(:) !! Seasonal MA coefficients by seasonal lag.
      real(dp), allocatable :: regression(:) !! User-regressor coefficients in input-column order.
      real(dp), allocatable :: residuals(:) !! CSS innovations or standardized ML residuals.
      real(dp), allocatable :: innovation_variance(:) !! Relative one-step innovation variances.
      real(dp), allocatable :: series(:) !! Original observations retained for forecasting.
      real(dp), allocatable :: adjusted_series(:) !! Observations after removing regression effects.
      real(dp), allocatable :: working_series(:) !! Series after ordinary and seasonal differences.
      logical, allocatable :: fixed(:) !! Mask in `[ar, ma, sar, sma, mean?, xreg]` order.
      real(dp), allocatable :: coef_covariance(:, :) !! Parameter covariance in the same full order.
      real(dp) :: mean = 0.0_dp !! Mean of the stationary series; zero when differencing is used.
      real(dp) :: variance = 0.0_dp !! Conditional residual mean square.
      real(dp) :: log_likelihood = 0.0_dp !! Gaussian CSS log likelihood reported by R.
      integer :: difference_order = 0 !! Number of ordinary differences applied before fitting.
      integer :: seasonal_difference_order = 0 !! Number of fitted seasonal differences.
      integer :: period = 1 !! Seasonal period, or one for a nonseasonal model.
      integer :: n_conditioned = 0 !! Number of leading original observations conditioned upon.
      integer :: n_used = 0 !! Number of observations after ordinary differencing.
      integer :: iterations = 0 !! Number of completed Nelder-Mead iterations.
      integer :: status = 0 !! Zero on convergence, one invalid input, or two iteration limit.
      integer :: covariance_status = 0 !! Zero when coefficient covariance was computed successfully.
      logical :: include_mean = .false. !! Whether a stationary-series mean was estimated.
      character(len=8) :: method = "CSS" !! Estimation method, currently `CSS` or `ML`.
   end type arima_fit_t

   type :: arima_forecast_t
      !! Stores ARIMA point forecasts and normal prediction intervals.
      real(dp), allocatable :: fit(:) !! Point forecasts for each requested horizon.
      real(dp), allocatable :: standard_error(:) !! Innovation-based forecast standard errors.
      real(dp), allocatable :: lower(:) !! Lower endpoints of the prediction intervals.
      real(dp), allocatable :: upper(:) !! Upper endpoints of the prediction intervals.
      real(dp) :: confidence_level = 0.95_dp !! Central interval probability.
      integer :: status = 0 !! Zero on success or one for invalid input.
   end type arima_forecast_t

   type :: holt_winters_fit_t
      !! Stores a Holt-Winters filter with explicit smoothing parameters and initial states.
      real(dp), allocatable :: fitted(:) !! One-step predictions after the initial period.
      real(dp), allocatable :: residuals(:) !! Observations minus one-step predictions.
      real(dp), allocatable :: level(:) !! Level state used for each one-step prediction.
      real(dp), allocatable :: trend(:) !! Trend state used for each one-step prediction.
      real(dp), allocatable :: seasonal(:) !! Seasonal state used for each one-step prediction.
      real(dp), allocatable :: final_season(:) !! Final seasonal states in forecast order.
      real(dp), allocatable :: initial_season(:) !! Seasonal states supplied to the filter.
      real(dp) :: final_level = 0.0_dp !! Level state after processing the final observation.
      real(dp) :: final_trend = 0.0_dp !! Trend state after processing the final observation.
      real(dp) :: initial_level = 0.0_dp !! Level state supplied to the filter.
      real(dp) :: initial_trend = 0.0_dp !! Trend state supplied to the filter.
      real(dp) :: sse = 0.0_dp !! Sum of squared one-step prediction errors.
      real(dp) :: alpha = 0.0_dp !! Level smoothing parameter.
      real(dp) :: beta = 0.0_dp !! Trend smoothing parameter.
      real(dp) :: gamma = 0.0_dp !! Seasonal smoothing parameter.
      integer :: period = 0 !! Number of observations in one seasonal cycle.
      integer :: iterations = 0 !! Optimization iterations, or zero for fixed parameters.
      integer :: status = 0 !! Zero on success; positive values identify invalid input or nonconvergence.
      logical :: has_season = .true. !! Whether the model contains a seasonal component.
      logical :: has_trend = .true. !! Whether the model contains a trend component.
      logical :: multiplicative = .false. !! Whether the seasonal recurrence is multiplicative.
      logical :: optimized = .false. !! Whether smoothing parameters were optimized.
   end type holt_winters_fit_t

   type :: holt_winters_forecast_t
      !! Stores Holt-Winters point forecasts and normal prediction intervals.
      real(dp), allocatable :: fit(:) !! Point forecasts for each requested horizon.
      real(dp), allocatable :: upper(:) !! Upper endpoints of the prediction intervals.
      real(dp), allocatable :: lower(:) !! Lower endpoints of the prediction intervals.
      real(dp), allocatable :: standard_error(:) !! Forecast standard errors.
      real(dp) :: confidence_level = 0.95_dp !! Central interval probability.
      integer :: status = 0 !! Zero on success or positive for invalid input.
   end type holt_winters_forecast_t

   type :: seasonal_decomposition_t
      !! Stores a classical moving-average seasonal decomposition.
      real(dp), allocatable :: trend(:) !! Moving-average trend; endpoints are quiet NaNs.
      real(dp), allocatable :: seasonal(:) !! Repeated normalized seasonal figure.
      real(dp), allocatable :: random(:) !! Remainder after removing trend and seasonality.
      real(dp), allocatable :: figure(:) !! Normalized seasonal pattern with shape `(period)`.
      integer :: period = 0 !! Number of observations in one seasonal cycle.
      integer :: status = 0 !! Zero on success or positive for invalid input.
      logical :: multiplicative = .false. !! Whether components combine multiplicatively.
   end type seasonal_decomposition_t

   type :: stl_result_t
      !! Stores a seasonal-trend decomposition based on LOESS.
      real(dp), allocatable :: seasonal(:) !! Locally smoothed seasonal component.
      real(dp), allocatable :: trend(:) !! Locally smoothed trend component.
      real(dp), allocatable :: remainder(:) !! Observations minus seasonal and trend components.
      real(dp), allocatable :: weights(:) !! Final robustness weights.
      integer :: window(3) = 0 !! Seasonal, trend, and low-pass smoothing windows.
      integer :: degree(3) = 0 !! Seasonal, trend, and low-pass local polynomial degrees.
      integer :: jump(3) = 0 !! Seasonal, trend, and low-pass smoothing jumps.
      integer :: inner = 0 !! Number of inner iterations per robustness pass.
      integer :: outer = 0 !! Number of robustness iterations.
      integer :: period = 0 !! Number of observations in one seasonal cycle.
      integer :: status = 0 !! Zero on success or positive for invalid input.
      logical :: periodic = .false. !! Whether seasonal values are constant within each phase.
   end type stl_result_t

   type :: spectrum_result_t
      !! Stores a univariate periodogram estimate and its Fourier frequencies.
      real(dp), allocatable :: frequency(:) !! Positive Fourier frequencies in cycles per unit time.
      real(dp), allocatable :: spectrum(:) !! One-sided spectral estimates at `frequency`.
      real(dp), allocatable :: kernel(:) !! Full smoothing kernel, or `[1]` when unsmoothed.
      real(dp) :: degrees_freedom = 0.0_dp !! Approximate spectral degrees of freedom.
      real(dp) :: bandwidth = 0.0_dp !! Equivalent spectral bandwidth in cycles per unit time.
      real(dp) :: taper = 0.0_dp !! Fraction tapered at each end of the input series.
      integer :: n_used = 0 !! Padded Fourier-transform length.
      integer :: n_original = 0 !! Number of input observations.
      integer :: order = 0 !! Fitted AR order for an autoregressive spectrum.
      integer :: status = 0 !! Zero on success or a shared-core invalid-input status.
      logical :: demeaned = .false. !! Whether a constant mean was removed explicitly.
      logical :: detrended = .false. !! Whether a linear trend and its intercept were removed.
      logical :: smoothed = .false. !! Whether a nontrivial spectral kernel was applied.
   end type spectrum_result_t

   type :: multivariate_spectrum_result_t
      !! Stores multivariate spectra and pairwise frequency-domain relationships.
      real(dp), allocatable :: frequency(:) !! Positive Fourier frequencies in cycles per unit time.
      real(dp), allocatable :: spectrum(:, :) !! Autospectra with shape `(frequency, series)`.
      complex(dp), allocatable :: cross_spectrum(:, :, :) !! Cross-spectra by frequency and series pair.
      real(dp), allocatable :: coherence(:, :, :) !! Squared coherence by frequency and series pair.
      real(dp), allocatable :: phase(:, :, :) !! Cross-spectrum phases in radians.
      real(dp), allocatable :: kernel(:) !! Full smoothing kernel, or `[1]` when unsmoothed.
      real(dp) :: degrees_freedom = 0.0_dp !! Approximate spectral degrees of freedom.
      real(dp) :: bandwidth = 0.0_dp !! Equivalent spectral bandwidth in cycles per unit time.
      real(dp) :: taper = 0.0_dp !! Fraction tapered at each end of every input series.
      integer :: n_used = 0 !! Padded Fourier-transform length.
      integer :: n_original = 0 !! Number of input observations per series.
      integer :: status = 0 !! Zero on success or a shared-core invalid-input status.
      logical :: demeaned = .false. !! Whether each series mean was removed explicitly.
      logical :: detrended = .false. !! Whether each series was linearly detrended.
      logical :: smoothed = .false. !! Whether a nontrivial spectral kernel was applied.
   end type multivariate_spectrum_result_t

   type :: smooth_xy_t
      !! Stores paired abscissae and fitted values returned by a smoother.
      real(dp), allocatable :: x(:) !! Sorted evaluation locations.
      real(dp), allocatable :: y(:) !! Smoothed values corresponding to `x`.
   end type smooth_xy_t

   type :: loess_fit_t
      !! Stores a univariate direct-surface local-polynomial regression fit.
      real(dp), allocatable :: x(:) !! Training predictor values in input order.
      real(dp), allocatable :: y(:) !! Training responses in input order.
      real(dp), allocatable :: weights(:) !! Nonnegative prior observation weights.
      real(dp), allocatable :: robustness(:) !! Final Tukey-bisquare robustness weights.
      real(dp), allocatable :: fitted_values(:) !! Fitted responses in input order.
      real(dp), allocatable :: residuals(:) !! Training responses minus fitted values.
      real(dp) :: span = 0.75_dp !! Fraction of observations defining each neighborhood.
      integer :: degree = 2 !! Local-polynomial degree from zero through two.
      integer :: iterations = 1 !! Number of local-regression fitting passes performed.
      integer :: status = 0 !! Zero on success or a shared invalid-input status.
      character(len=9) :: family = "gaussian" !! `gaussian` or robust `symmetric` fitting.
   end type loess_fit_t

   type :: smooth_spline_fit_t
      !! Stores an all-knots univariate cubic smoothing-spline fit.
      real(dp), allocatable :: x(:) !! Sorted distinct training predictor values.
      real(dp), allocatable :: y(:) !! Training responses in sorted predictor order.
      real(dp), allocatable :: weights(:) !! Normalized nonnegative observation weights.
      real(dp), allocatable :: fitted_values(:) !! Smoothed responses at the training points.
      real(dp), allocatable :: residuals(:) !! Training responses minus fitted values.
      real(dp), allocatable :: leverage(:) !! Diagonal elements of the smoothing matrix.
      real(dp), allocatable :: knots(:) !! Normalized clamped cubic B-spline knots.
      real(dp), allocatable :: coefficients(:) !! Cubic B-spline coefficients.
      real(dp) :: minimum_x = 0.0_dp !! Lower endpoint used to normalize predictors.
      real(dp) :: x_range = 0.0_dp !! Positive training-predictor range.
      real(dp) :: degrees_freedom = 0.0_dp !! Trace of the smoothing matrix.
      real(dp) :: spar = 0.0_dp !! R-compatible smoothing parameter.
      real(dp) :: lambda = 0.0_dp !! Roughness-penalty multiplier.
      real(dp) :: ratio = 0.0_dp !! Design-to-roughness trace ratio used to map `spar`.
      real(dp) :: gcv = 0.0_dp !! Generalized cross-validation criterion at the fitted solution.
      real(dp) :: cv = 0.0_dp !! Observation-level leave-one-out cross-validation score.
      real(dp) :: criterion = 0.0_dp !! Criterion corresponding to the requested selection mode.
      real(dp) :: tolerance = 0.0_dp !! Predictor-grouping tolerance used by the fit.
      integer :: n_original = 0 !! Input count before duplicate aggregation.
      integer :: n_knots = 0 !! Number of distinct predictor knots used by the spline basis.
      integer :: status = 0 !! Zero on success or a shared numerical status.
      logical :: automatic = .false. !! True when a CV criterion selected `spar`.
      logical :: ordinary_cv = .false. !! True when ordinary CV replaced GCV as the criterion.
   end type smooth_spline_fit_t

   type :: nls_fit_t
      !! Stores a formula-free nonlinear least-squares fit.
      real(dp), allocatable :: coefficients(:) !! Estimated model parameters.
      real(dp), allocatable :: fitted_values(:) !! Model values at the solution.
      real(dp), allocatable :: residuals(:) !! Observations minus fitted values.
      real(dp), allocatable :: covariance(:, :) !! Estimated parameter covariance matrix.
      real(dp), allocatable :: standard_errors(:) !! Parameter standard errors.
      real(dp), allocatable :: lower_bounds(:) !! Applied lower parameter bounds.
      real(dp), allocatable :: upper_bounds(:) !! Applied upper parameter bounds.
      real(dp) :: rss = 0.0_dp !! Residual sum of squares.
      real(dp) :: residual_standard_error = 0.0_dp !! Square root of RSS divided by residual df.
      integer :: degrees_freedom = 0 !! Residual degrees of freedom.
      integer :: iterations = 0 !! Accepted Gauss-Newton iterations.
      integer :: rank = 0 !! Numerical rank of the final model Jacobian.
      integer :: status = 0 !! Zero on convergence or a numerical status code.
      logical :: converged = .false. !! True when the convergence tolerance was met.
      logical :: analytic_jacobian = .false. !! True when the caller supplied a Jacobian.
      logical :: bounded = .false. !! True when at least one parameter bound was supplied.
   end type nls_fit_t

   type :: kmeans_result_t
      !! Stores a deterministic Lloyd k-means fit initialized by explicit centers.
      real(dp), allocatable :: centers(:, :) !! Final cluster centers with shape `(k, p)`.
      integer, allocatable :: cluster(:) !! One-based cluster assignment for each observation.
      integer, allocatable :: size(:) !! Number of observations assigned to each cluster.
      real(dp), allocatable :: withinss(:) !! Within-cluster sums of squared distances.
      real(dp) :: totss = 0.0_dp !! Total sum of squares around the overall column means.
      real(dp) :: tot_withinss = 0.0_dp !! Sum of all within-cluster sums of squares.
      real(dp) :: betweenss = 0.0_dp !! Difference between total and within-cluster sums of squares.
      integer :: iter = 0 !! Number of Lloyd assignment/update iterations performed.
      integer :: status = 0 !! Zero on convergence, one invalid input, two iteration limit, or three empty cluster.
   end type kmeans_result_t

   type :: prcomp_fit_t
      !! Stores a principal-components fit corresponding to R's `prcomp` result.
      real(dp), allocatable :: sdev(:) !! Component standard deviations in descending order.
      real(dp), allocatable :: rotation(:, :) !! Variable loadings with shape `(p, rank)`.
      real(dp), allocatable :: x(:, :) !! Observation scores with shape `(n, rank)`.
      real(dp), allocatable :: center(:) !! Column centers applied before decomposition.
      real(dp), allocatable :: scale(:) !! Column scales applied before decomposition.
      integer :: status = 0 !! Zero on success, or a standardization or LAPACK status.
   end type prcomp_fit_t

   type :: scale_result_t
      !! Stores a standardized matrix and the column transformations applied.
      real(dp), allocatable :: values(:, :) !! Centered and/or scaled data matrix.
      real(dp), allocatable :: center(:) !! Applied column centers.
      real(dp), allocatable :: scale(:) !! Applied column divisors.
      integer :: status = 0 !! Zero on success or a shared-core status code.
   end type scale_result_t

   type :: fisher_test_result_t
      !! Stores a two-sided Fisher exact test for a two-by-two table.
      real(dp) :: p_value = 1.0_dp !! Two-sided conditional exact p-value.
      real(dp) :: estimate = 0.0_dp !! Conditional maximum-likelihood odds-ratio estimate.
      integer :: method = 1 !! Method identifier reserved for future extensions.
   end type fisher_test_result_t

   type :: kruskal_test_result_t
      !! Stores a Kruskal-Wallis rank-sum test.
      real(dp) :: statistic = 0.0_dp !! Tie-corrected Kruskal-Wallis chi-square statistic.
      integer :: parameter = 0 !! Number of nonempty groups minus one.
      real(dp) :: p_value = 1.0_dp !! Upper-tail chi-square approximation.
   end type kruskal_test_result_t

   type :: ks_test_result_t
      !! Stores an asymptotic one- or two-sample Kolmogorov-Smirnov test.
      real(dp) :: statistic = 0.0_dp !! Maximum absolute empirical-distribution difference.
      real(dp) :: p_value = 1.0_dp !! Two-sided asymptotic p-value.
      integer :: n = 0 !! One-sample size or combined two-sample size.
      integer :: method = 1 !! One for a normal reference and two for a two-sample test.
   end type ks_test_result_t

   type :: wilcox_test_result_t
      !! Stores an asymptotic Wilcoxon rank-sum or signed-rank test.
      real(dp) :: statistic = 0.0_dp !! Rank-sum W statistic.
      real(dp) :: p_value = 1.0_dp !! Two-sided normal-approximation p-value.
      integer :: method = 1 !! One for rank-sum and two for paired signed-rank.
   end type wilcox_test_result_t

   type :: chisq_test_result_t
      !! Stores a Pearson chi-square goodness-of-fit or independence test.
      real(dp) :: statistic = 0.0_dp !! Pearson chi-square statistic, with Yates correction when selected.
      integer :: parameter = 0 !! Degrees of freedom.
      real(dp) :: p_value = 1.0_dp !! Upper-tail chi-square p-value.
      integer :: method = 1 !! One=goodness of fit, two=independence, three=Bartlett, four=Friedman.
   end type chisq_test_result_t

   type :: cor_test_result_t
      !! Stores a Pearson or asymptotic Spearman correlation test.
      real(dp) :: statistic = 0.0_dp !! Pearson t statistic or Spearman S statistic.
      integer :: parameter = 0 !! Pearson degrees of freedom, also retained for asymptotic Spearman.
      real(dp) :: p_value = 1.0_dp !! Two-sided p-value.
      real(dp) :: estimate = 0.0_dp !! Pearson correlation or Spearman rank correlation.
      real(dp) :: conf_low = 0.0_dp !! Lower Pearson confidence limit, when available.
      real(dp) :: conf_high = 0.0_dp !! Upper Pearson confidence limit, when available.
      integer :: method = 1 !! One for Pearson and two for asymptotic Spearman.
   end type cor_test_result_t

   type :: prop_test_result_t
      !! Stores a one- or multiple-sample Pearson proportion test.
      real(dp) :: statistic = 0.0_dp !! Pearson chi-square statistic.
      integer :: parameter = 0 !! Degrees of freedom.
      real(dp) :: p_value = 1.0_dp !! Two-sided p-value.
      real(dp) :: estimate = 0.0_dp !! First estimated sample proportion.
      real(dp) :: estimate2 = 0.0_dp !! Second estimated sample proportion, when present.
      real(dp), allocatable :: estimates(:) !! All estimated sample proportions.
      real(dp) :: null_value = 0.0_dp !! Scalar null proportion or pooled estimate.
      real(dp), allocatable :: null_values(:) !! Null proportions when explicitly supplied.
      real(dp) :: conf_low = 0.0_dp !! Lower 95 percent confidence limit when defined.
      real(dp) :: conf_high = 0.0_dp !! Upper 95 percent confidence limit when defined.
      integer :: method = 1 !! One for a one-sample test and two for a multiple-sample test.
   end type prop_test_result_t

   type :: t_test_result_t
      !! Stores numerical results from a one-sample, two-sample, or paired t test.
      real(dp) :: statistic = 0.0_dp !! Student t statistic.
      real(dp) :: parameter = 0.0_dp !! Degrees of freedom, possibly nonintegral for Welch's test.
      real(dp) :: p_value = 1.0_dp !! Two-sided p-value.
      real(dp) :: estimate = 0.0_dp !! Sample mean or first sample mean.
      real(dp) :: estimate2 = 0.0_dp !! Second sample mean for an independent two-sample test.
      real(dp) :: null_value = 0.0_dp !! Null-hypothesis value for the tested mean or mean difference.
      real(dp) :: conf_low = 0.0_dp !! Lower endpoint of the 95 percent confidence interval.
      real(dp) :: conf_high = 0.0_dp !! Upper endpoint of the 95 percent confidence interval.
      real(dp) :: stderr = 0.0_dp !! Standard error of the estimated mean or mean difference.
      integer :: method = 1 !! One=one-sample, two=Welch, three=pooled, and four=paired.
   end type t_test_result_t

   type :: glm_fit_t
      !! Stores numerical results and fitting data for a supported generalized linear model.
      real(dp), allocatable :: coef(:) !! Estimated coefficients, including the intercept.
      real(dp), allocatable :: se(:) !! Asymptotic coefficient standard errors.
      real(dp), allocatable :: z_value(:) !! Coefficient Wald z statistics.
      real(dp), allocatable :: p_value(:) !! Two-sided normal-approximation coefficient p-values.
      real(dp), allocatable :: fitted(:) !! Fitted response means with size `n`.
      real(dp), allocatable :: resid(:) !! Response residuals with size `n`.
      real(dp), allocatable :: y(:) !! Response values used to fit the model.
      real(dp), allocatable :: xpred(:, :) !! Predictor matrix without an intercept column.
      real(dp), allocatable :: offset(:) !! Additive linear-predictor offset used during fitting.
      real(dp), allocatable :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), allocatable :: working_weights(:) !! Final IRLS working weights with size `n`.
      real(dp), allocatable :: covariance(:, :) !! Estimated coefficient covariance matrix.
      logical, allocatable :: aliased(:) !! True for coefficients excluded by rank detection.
      real(dp) :: dispersion = 1.0_dp !! Estimated or fixed family dispersion.
      real(dp) :: deviance = 0.0_dp !! Residual deviance at the fitted means.
      integer :: df = 0 !! Residual degrees of freedom.
      integer :: rank = 0 !! Numerical rank of the weighted design matrix.
      integer :: convergence = 1 !! Zero after successful convergence and one otherwise.
      integer :: iter = 0 !! Number of Newton iterations performed.
      integer :: family = 0 !! Family identifier defined by `r_stats_glm`.
      integer :: link = 0 !! Link identifier defined by `r_stats_glm`.
      logical :: has_intercept = .true. !! Whether the fitted design includes an intercept.
      logical :: weighted = .false. !! Whether caller-supplied prior weights were used.
   end type glm_fit_t

   type :: ecdf_t
      !! Stores the knots and cumulative probabilities of an empirical distribution function.
      real(dp), allocatable :: values(:) !! Sorted distinct sample values.
      real(dp), allocatable :: probabilities(:) !! Cumulative probabilities at the sample values.
      integer :: n = 0 !! Number of observations represented, including duplicates.
   end type ecdf_t

   type :: histogram_t
      !! Stores a numeric histogram without plotting state.
      real(dp), allocatable :: breaks(:) !! Strictly increasing bin boundaries.
      real(dp), allocatable :: mids(:) !! Bin midpoints.
      real(dp), allocatable :: density(:) !! Bin densities integrating to one.
      integer, allocatable :: counts(:) !! Observation counts in each bin.
      logical :: equidistant = .false. !! True when all bin widths are equal within rounding error.
   end type histogram_t

   type :: density_result_t
      !! Stores a univariate kernel-density estimate evaluated on a regular grid.
      real(dp), allocatable :: x(:) !! Increasing evaluation grid.
      real(dp), allocatable :: y(:) !! Estimated density at each grid point.
      real(dp) :: bandwidth = 0.0_dp !! Applied kernel bandwidth.
      integer :: n = 0 !! Number of observations represented.
      integer :: kernel = 0 !! Kernel identifier defined by `r_stats_empirical`.
   end type density_result_t

   type :: interpolation_t
      !! Stores sorted knots and controls for reusable linear or constant interpolation.
      real(dp), allocatable :: x(:) !! Strictly increasing interpolation knots.
      real(dp), allocatable :: y(:) !! Responses at the knots after resolving ties.
      real(dp) :: constant_mix = 0.0_dp !! Constant interpolation mix corresponding to R's `f`.
      integer :: method = 1 !! One for linear and two for constant interpolation.
      integer :: left_rule = 1 !! One for NaN and two for constant left extrapolation.
      integer :: right_rule = 1 !! One for NaN and two for constant right extrapolation.
   end type interpolation_t

   type :: cubic_spline_t
      !! Stores piecewise-cubic coefficients for a reusable interpolation spline.
      real(dp), allocatable :: x(:) !! Strictly increasing spline knots.
      real(dp), allocatable :: a(:) !! Constant coefficient on each knot interval.
      real(dp), allocatable :: b(:) !! Linear coefficient on each knot interval.
      real(dp), allocatable :: c(:) !! Quadratic coefficient on each knot interval.
      real(dp), allocatable :: d(:) !! Cubic coefficient on each knot interval.
      real(dp) :: left_slope = 0.0_dp !! Linear extrapolation slope for a natural spline.
      real(dp) :: right_slope = 0.0_dp !! Linear extrapolation slope for a natural spline.
      integer :: method = 1 !! One for FMM and two for natural boundary conditions.
   end type cubic_spline_t

   type :: isoreg_t
      !! Stores an increasing isotonic least-squares fit.
      real(dp), allocatable :: x(:) !! Increasing predictor values.
      real(dp), allocatable :: y(:) !! Responses in sorted predictor order.
      real(dp), allocatable :: fitted(:) !! Monotone fitted responses.
      real(dp), allocatable :: weights(:) !! Nonnegative fitting weights.
      integer, allocatable :: knots(:) !! Final indices of constant fitted blocks.
   end type isoreg_t

   type :: lm_fit_t
      !! Stores the numerical results and design data for a fitted linear model.
      real(dp), allocatable :: coef(:) !! Estimated coefficients, including the intercept when requested.
      real(dp), allocatable :: fitted(:) !! Fitted response values with size `n`.
      real(dp), allocatable :: resid(:) !! Response residuals with size `n`.
      real(dp), allocatable :: cov_unscaled(:, :) !! Unscaled coefficient covariance, `(X**T X)**-1`.
      real(dp), allocatable :: y(:) !! Response values used to fit the model.
      real(dp), allocatable :: xpred(:, :) !! Predictor matrix without an intercept column.
      real(dp), allocatable :: weights(:) !! Nonnegative prior weights with size `n`.
      logical, allocatable :: aliased(:) !! True for coefficients excluded by rank detection.
      real(dp) :: sigma = 0.0_dp !! Residual standard error.
      real(dp) :: r_squared = 0.0_dp !! Coefficient of determination.
      real(dp) :: adj_r_squared = 0.0_dp !! Degrees-of-freedom-adjusted coefficient of determination.
      integer :: df = 0 !! Residual degrees of freedom.
      integer :: rank = 0 !! Numerical rank of the weighted design matrix.
      logical :: has_intercept = .true. !! Whether the fitted design includes an intercept.
      logical :: weighted = .false. !! Whether the fit used caller-supplied prior weights.
   end type lm_fit_t

   type :: one_way_anova_t
      !! Stores a formula-free one-way analysis of variance.
      real(dp), allocatable :: group_means(:) !! Mean response for each integer group.
      integer, allocatable :: group_counts(:) !! Observation count for each group.
      real(dp) :: ss_between = 0.0_dp !! Between-group sum of squares.
      real(dp) :: ss_within = 0.0_dp !! Within-group sum of squares.
      real(dp) :: ms_between = 0.0_dp !! Between-group mean square.
      real(dp) :: ms_within = 0.0_dp !! Within-group mean square.
      real(dp) :: f_statistic = 0.0_dp !! Omnibus F statistic.
      real(dp) :: p_value = 1.0_dp !! Upper-tail F probability.
      integer :: df_between = 0 !! Between-group degrees of freedom.
      integer :: df_within = 0 !! Within-group degrees of freedom.
      integer :: status = 0 !! Zero on success or nonzero for invalid input.
   end type one_way_anova_t

   type :: anova_comparison_t
      !! Stores an extra-sum-of-squares comparison of nested linear models.
      real(dp) :: sum_of_squares = 0.0_dp !! Reduction in residual sum of squares.
      real(dp) :: mean_square = 0.0_dp !! Reduction divided by numerator degrees of freedom.
      real(dp) :: f_statistic = 0.0_dp !! Partial F statistic.
      real(dp) :: p_value = 1.0_dp !! Upper-tail F probability.
      integer :: df_numerator = 0 !! Difference in fitted model ranks.
      integer :: df_denominator = 0 !! Residual degrees of freedom of the larger model.
      integer :: status = 0 !! Zero on success or nonzero for invalid nesting.
   end type anova_comparison_t

   type :: mlm_fit_t
      !! Stores a full-rank multivariate linear-model fit.
      real(dp), allocatable :: coefficients(:, :) !! Coefficients with shape `(k,m)`.
      real(dp), allocatable :: fitted(:, :) !! Fitted responses with shape `(n,m)`.
      real(dp), allocatable :: residuals(:, :) !! Response residuals with shape `(n,m)`.
      real(dp), allocatable :: residual_sscp(:, :) !! Residual SSCP matrix with shape `(m,m)`.
      real(dp), allocatable :: y(:, :) !! Response matrix used for fitting.
      integer :: rank = 0 !! Design-matrix rank.
      integer :: df_residual = 0 !! Residual degrees of freedom.
      integer :: status = 0 !! Zero on success or nonzero on invalid input or failure.
      logical :: has_intercept = .true. !! Whether the model includes an intercept.
   end type mlm_fit_t

   type :: manova_result_t
      !! Stores four standard tests comparing nested multivariate linear models.
      real(dp) :: statistic(4) = 0.0_dp !! Pillai, Wilks, Hotelling-Lawley, and Roy criteria.
      real(dp) :: approximate_f(4) = 0.0_dp !! Corresponding approximate F statistics.
      real(dp) :: numerator_df(4) = 0.0_dp !! Approximate numerator degrees of freedom.
      real(dp) :: denominator_df(4) = 0.0_dp !! Approximate denominator degrees of freedom.
      real(dp) :: p_value(4) = 1.0_dp !! Upper-tail F probabilities.
      real(dp), allocatable :: eigenvalues(:) !! Eigenvalues of `E**-1 H`.
      integer :: hypothesis_df = 0 !! Difference in model ranks.
      integer :: status = 0 !! Zero on success or nonzero for invalid models.
   end type manova_result_t

   type :: cancor_result_t
      !! Stores a canonical-correlation analysis.
      real(dp), allocatable :: correlation(:) !! Canonical correlations in decreasing order.
      real(dp), allocatable :: x_coefficients(:, :) !! Canonical coefficients for `x`.
      real(dp), allocatable :: y_coefficients(:, :) !! Canonical coefficients for `y`.
      real(dp), allocatable :: x_center(:) !! Centers subtracted from `x`.
      real(dp), allocatable :: y_center(:) !! Centers subtracted from `y`.
      integer :: status = 0 !! Zero on success or nonzero for invalid or singular input.
   end type cancor_result_t

end module r_stats_types
