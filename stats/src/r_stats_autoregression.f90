! SPDX-License-Identifier: MIT
! SPDX-FileComment: Autoregressive fitting and spectra corresponding to selected R stats functions.
module r_stats_autoregression
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_linalg, only: solve_system
   use r_missing, only: r_is_finite
   use r_optional, only: optval
   use r_status, only: r_invalid_input, r_ok
   use r_time_series, only: r_arma_autocorrelation, r_autocovariance, r_durbin_levinson
   use r_stats_types, only: ar_fit_t, spectrum_result_t
   implicit none
   private

   public :: ar, ar_burg, ar_mle, ar_ols, ar_yw, spec_ar

contains

   pure function ar_yw(x, order, order_max, demean) result(fit)
      !! Fits a univariate autoregression by the Yule-Walker method.
      real(dp), intent(in) :: x(:) !! Finite univariate observations.
      integer, intent(in), optional :: order !! Fixed AR order; omission selects an order by AIC.
      integer, intent(in), optional :: order_max !! Largest AIC candidate order.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to true.
      type(ar_fit_t) :: fit
      real(dp), allocatable :: ar_table(:, :), autocovariance(:), partial(:), variances(:)
      real(dp) :: mean_value, minimum_aic
      integer :: i, j, maximum_order, selected_order
      logical :: remove_mean

      if (size(x) < 2 .or. any(.not. r_is_finite(x))) then
         fit%status = r_invalid_input
         return
      end if
      remove_mean = optval(demean, .true.)
      if (present(order)) then
         maximum_order = order
      else if (present(order_max)) then
         maximum_order = order_max
      else
         maximum_order = min(size(x) - 1, int(floor(10.0_dp*log10(real(size(x), dp)))))
      end if
      if (maximum_order < 0 .or. maximum_order >= size(x)) then
         fit%status = r_invalid_input
         return
      end if

      call r_autocovariance(x, autocovariance, lag_max=maximum_order, &
                            demean=remove_mean, status=fit%status)
      if (fit%status /= r_ok) return
      call r_durbin_levinson(autocovariance, partial, fit%status, ar_table, variances)
      if (fit%status /= r_ok) return
      allocate (fit%aic(maximum_order + 1))
      allocate (fit%aic_order(maximum_order + 1))
      do i = 0, maximum_order
         fit%aic_order(i + 1) = i
         fit%aic(i + 1) = real(size(x), dp)*log(variances(i)) + 2.0_dp*real(i, dp)
      end do
      minimum_aic = minval(fit%aic)
      fit%aic = fit%aic - minimum_aic
      if (present(order)) then
         selected_order = order
      else
         selected_order = minloc(fit%aic, dim=1) - 1
      end if
      allocate (fit%coefficients(selected_order))
      if (selected_order > 0) fit%coefficients = ar_table(selected_order, 1:selected_order)
      mean_value = 0.0_dp
      if (remove_mean) mean_value = sum(x)/real(size(x), dp)
      fit%mean = mean_value
      fit%order = selected_order
      fit%variance = variances(selected_order)*real(size(x), dp)/ &
                     real(size(x) - selected_order - merge(1, 0, remove_mean), dp)
      call calculate_residuals(x, fit%coefficients, mean_value, 0.0_dp, fit%residuals)
      fit%n_used = size(x)
      fit%demeaned = remove_mean
   end function ar_yw

   pure function ar_burg(x, order, order_max, demean, variance_method) result(fit)
      !! Fits a univariate autoregression with Burg's forward-backward recursion.
      real(dp), intent(in) :: x(:) !! Finite univariate observations.
      integer, intent(in), optional :: order !! Fixed AR order; omission selects an order by AIC.
      integer, intent(in), optional :: order_max !! Largest AIC candidate order.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to true.
      integer, intent(in), optional :: variance_method !! One for recursive, two for residual.
      type(ar_fit_t) :: fit
      real(dp), allocatable :: coefficients(:, :), forward(:), backward(:), old_forward(:)
      real(dp), allocatable :: variance_one(:), variance_two(:), y(:)
      real(dp) :: denominator, mean_value, reflection, sum_squares
      integer :: i, j, maximum_order, method_number, selected_order
      logical :: remove_mean

      if (size(x) < 2 .or. any(.not. r_is_finite(x))) then
         fit%status = r_invalid_input
         return
      end if
      remove_mean = optval(demean, .true.)
      method_number = optval(variance_method, 1)
      if (present(order)) then
         maximum_order = order
      else if (present(order_max)) then
         maximum_order = order_max
      else
         maximum_order = min(size(x) - 1, int(floor(10.0_dp*log10(real(size(x), dp)))))
      end if
      if (maximum_order < 0 .or. maximum_order >= size(x) .or. &
          (method_number /= 1 .and. method_number /= 2)) then
         fit%status = r_invalid_input
         return
      end if
      mean_value = 0.0_dp
      if (remove_mean) mean_value = sum(x)/real(size(x), dp)
      y = x - mean_value
      sum_squares = sum(y**2)
      if (sum_squares <= tiny(1.0_dp)) then
         fit%status = r_invalid_input
         return
      end if
      allocate (forward(size(x)), backward(size(x)), old_forward(size(x)))
      allocate (coefficients(maximum_order, maximum_order), source=0.0_dp)
      allocate (variance_one(0:maximum_order), variance_two(0:maximum_order))
      do i = 1, size(x)
         forward(i) = y(size(x) + 1 - i)
      end do
      backward = forward
      variance_one(0) = sum_squares/real(size(x), dp)
      variance_two(0) = variance_one(0)
      do i = 1, maximum_order
         sum_squares = 0.0_dp
         denominator = 0.0_dp
         do j = i + 1, size(x)
            sum_squares = sum_squares + backward(j)*forward(j - 1)
            denominator = denominator + backward(j)**2 + forward(j - 1)**2
         end do
         if (denominator <= tiny(1.0_dp)) then
            fit%status = r_invalid_input
            return
         end if
         reflection = 2.0_dp*sum_squares/denominator
         coefficients(i, i) = reflection
         do j = 1, i - 1
            coefficients(i, j) = coefficients(i - 1, j) - &
                                  reflection*coefficients(i - 1, i - j)
         end do
         old_forward = forward
         do j = i + 1, size(x)
            forward(j) = old_forward(j - 1) - reflection*backward(j)
            backward(j) = backward(j) - reflection*old_forward(j - 1)
         end do
         variance_one(i) = variance_one(i - 1)*(1.0_dp - reflection**2)
         sum_squares = 0.0_dp
         do j = i + 1, size(x)
            sum_squares = sum_squares + backward(j)**2 + forward(j)**2
         end do
         variance_two(i) = sum_squares/(2.0_dp*real(size(x) - i, dp))
      end do
      allocate (fit%aic(maximum_order + 1), fit%aic_order(maximum_order + 1))
      do i = 0, maximum_order
         fit%aic_order(i + 1) = i
         fit%aic(i + 1) = real(size(x), dp)*log(merge(variance_one(i), variance_two(i), &
                                                      method_number == 1)) + 2.0_dp*real(i, dp)
      end do
      fit%aic = fit%aic - minval(fit%aic)
      if (present(order)) then
         selected_order = order
      else
         selected_order = minloc(fit%aic, dim=1) - 1
      end if
      allocate (fit%coefficients(selected_order))
      if (selected_order > 0) fit%coefficients = coefficients(selected_order, 1:selected_order)
      fit%variance = merge(variance_one(selected_order), variance_two(selected_order), &
                           method_number == 1)
      fit%mean = mean_value
      fit%order = selected_order
      fit%n_used = size(x)
      fit%demeaned = remove_mean
      fit%method = "burg"
      call calculate_residuals(x, fit%coefficients, mean_value, 0.0_dp, fit%residuals)
   end function ar_burg

   pure function ar_ols(x, order, order_max, demean, intercept) result(fit)
      !! Fits a univariate autoregression by unconstrained least squares.
      real(dp), intent(in) :: x(:) !! Finite univariate observations.
      integer, intent(in), optional :: order !! Fixed AR order; omission selects an order by AIC.
      integer, intent(in), optional :: order_max !! Largest AIC candidate order.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to true.
      logical, intent(in), optional :: intercept !! Fit an intercept; defaults to `demean`.
      type(ar_fit_t) :: fit
      real(dp), allocatable :: beta(:), coefficients(:, :), design(:, :), response(:)
      real(dp), allocatable :: residuals(:), variances(:), y(:)
      real(dp) :: mean_value
      integer :: candidate_count, column, first_order, i, info, maximum_order, parameter_count
      integer :: selected_order
      logical :: include_intercept, remove_mean

      if (size(x) < 2 .or. any(.not. r_is_finite(x))) then
         fit%status = r_invalid_input
         return
      end if
      remove_mean = optval(demean, .true.)
      include_intercept = remove_mean
      if (present(intercept)) include_intercept = intercept
      if (present(order)) then
         maximum_order = order
         first_order = order
      else if (present(order_max)) then
         maximum_order = order_max
         first_order = 0
      else
         maximum_order = min(size(x) - 1, int(floor(10.0_dp*log10(real(size(x), dp)))))
         first_order = 0
      end if
      if (maximum_order < 0 .or. maximum_order >= size(x)) then
         fit%status = r_invalid_input
         return
      end if
      mean_value = 0.0_dp
      if (remove_mean) mean_value = sum(x)/real(size(x), dp)
      y = x - mean_value
      candidate_count = maximum_order - first_order + 1
      allocate (fit%aic(candidate_count), fit%aic_order(candidate_count))
      allocate (coefficients(maximum_order, maximum_order), source=0.0_dp)
      allocate (variances(0:maximum_order), source=huge(1.0_dp))
      do i = first_order, maximum_order
         parameter_count = i + merge(1, 0, include_intercept)
         allocate (response(size(x) - i))
         response = y(i + 1:)
         if (parameter_count > 0) then
            allocate (design(size(x) - i, parameter_count), beta(parameter_count))
            if (include_intercept) design(:, 1) = 1.0_dp
            do column = 1, i
               design(:, merge(1, 0, include_intercept) + column) = &
                  y(i + 1 - column:size(x) - column)
            end do
            call solve_system(matmul(transpose(design), design), &
                              matmul(transpose(design), response), beta, info)
            if (info /= 0) then
               fit%status = r_invalid_input
               return
            end if
            residuals = response - matmul(design, beta)
            if (include_intercept) fit%intercept = beta(1)
            if (i > 0) coefficients(i, 1:i) = beta(merge(1, 0, include_intercept) + 1:)
            deallocate (design, beta)
         else
            residuals = response
         end if
         variances(i) = sum(residuals**2)/real(size(x) - i, dp)
         fit%aic_order(i - first_order + 1) = i
         fit%aic(i - first_order + 1) = real(size(x), dp)*log(variances(i)) + &
                                        2.0_dp*real(parameter_count, dp)
         deallocate (response, residuals)
      end do
      fit%aic = fit%aic - minval(fit%aic)
      selected_order = fit%aic_order(minloc(fit%aic, dim=1))
      allocate (fit%coefficients(selected_order))
      if (selected_order > 0) fit%coefficients = coefficients(selected_order, 1:selected_order)
      fit%order = selected_order
      fit%variance = variances(selected_order)
      fit%mean = mean_value
      fit%n_used = size(x)
      fit%demeaned = remove_mean
      fit%method = "ols"
      call calculate_ols_intercept(y, selected_order, include_intercept, &
                                   fit%coefficients, fit%intercept)
      call calculate_residuals(x, fit%coefficients, mean_value, fit%intercept, fit%residuals)
   end function ar_ols

   pure function ar_mle(x, order, order_max, demean, delta) result(fit)
      !! Fits a stationary univariate Gaussian autoregression by exact maximum likelihood.
      real(dp), intent(in) :: x(:) !! Finite univariate observations.
      integer, intent(in), optional :: order !! Fixed AR order; omission selects an order by AIC.
      integer, intent(in), optional :: order_max !! Largest AIC candidate order.
      logical, intent(in), optional :: demean !! Estimate the stationary mean; defaults to true.
      real(dp), intent(in), optional :: delta !! Steady-state switch tolerance; defaults to zero.
      type(ar_fit_t) :: fit
      real(dp), allocatable :: all_coefficients(:, :), candidate_aic(:), candidate_mean(:)
      real(dp), allocatable :: candidate_variance(:), coefficients(:)
      real(dp) :: log_likelihood
      integer :: candidate_count, first_order, i, maximum_order, selected_index
      logical :: estimate_mean

      if (size(x) < 2 .or. any(.not. r_is_finite(x))) then
         fit%status = r_invalid_input
         return
      end if
      estimate_mean = optval(demean, .true.)
      if (optval(delta, 0.0_dp) < 0.0_dp) then
         fit%status = r_invalid_input
         return
      end if
      if (present(order)) then
         maximum_order = order
         first_order = order
      else if (present(order_max)) then
         maximum_order = order_max
         first_order = 0
      else
         maximum_order = min(size(x) - 1, 12, &
            int(floor(10.0_dp*log10(real(size(x), dp)))))
         first_order = 0
      end if
      if (maximum_order < 0 .or. maximum_order >= size(x)) then
         fit%status = r_invalid_input
         return
      end if
      candidate_count = maximum_order - first_order + 1
      allocate (candidate_aic(candidate_count), candidate_mean(candidate_count))
      allocate (candidate_variance(candidate_count))
      allocate (all_coefficients(maximum_order, candidate_count), source=0.0_dp)
      do i = first_order, maximum_order
         call fit_mle_order(x, i, estimate_mean, optval(delta, 0.0_dp), coefficients, &
                            candidate_mean(i - first_order + 1), &
                            candidate_variance(i - first_order + 1), log_likelihood, fit%status)
         if (fit%status /= r_ok) return
         if (i > 0) all_coefficients(1:i, i - first_order + 1) = coefficients
         candidate_aic(i - first_order + 1) = -2.0_dp*log_likelihood + &
            2.0_dp*real(i + merge(1, 0, estimate_mean) + 1, dp)
      end do
      selected_index = minloc(candidate_aic, dim=1)
      fit%order = selected_index + first_order - 1
      allocate (fit%coefficients(fit%order))
      if (fit%order > 0) fit%coefficients = all_coefficients(1:fit%order, selected_index)
      fit%mean = candidate_mean(selected_index)
      fit%variance = candidate_variance(selected_index)
      call stationary_ar_likelihood(x, fit%coefficients, estimate_mean, &
         optval(delta, 0.0_dp), fit%mean, fit%variance, fit%log_likelihood, fit%status)
      if (fit%status /= r_ok) return
      allocate (fit%aic(candidate_count), fit%aic_order(candidate_count))
      fit%aic = candidate_aic - minval(candidate_aic)
      do i = 1, candidate_count
         fit%aic_order(i) = first_order + i - 1
      end do
      fit%n_used = size(x)
      fit%demeaned = estimate_mean
      fit%method = "mle"
      call calculate_residuals(x, fit%coefficients, fit%mean, 0.0_dp, fit%residuals)
   end function ar_mle

   pure function ar(x, order, order_max, demean, method, intercept, variance_method, &
                    mle_delta) result(fit)
      !! Fits an autoregression using the requested supported R method.
      real(dp), intent(in) :: x(:) !! Finite univariate observations.
      integer, intent(in), optional :: order !! Fixed AR order; omission selects an order by AIC.
      integer, intent(in), optional :: order_max !! Largest AIC candidate order.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to true.
      character(len=*), intent(in), optional :: method !! `yule-walker`, `burg`, `ols`, or `mle`.
      logical, intent(in), optional :: intercept !! OLS intercept; defaults to `demean`.
      integer, intent(in), optional :: variance_method !! Burg variance method one or two.
      real(dp), intent(in), optional :: mle_delta !! MLE steady-state tolerance; defaults to zero.
      type(ar_fit_t) :: fit
      character(len=16) :: selected_method

      selected_method = "yule-walker"
      if (present(method)) selected_method = lowercase(method)
      select case (trim(selected_method))
      case ("yule-walker", "yw")
         fit = ar_yw(x, order, order_max, demean)
      case ("burg")
         fit = ar_burg(x, order, order_max, demean, variance_method)
      case ("ols")
         fit = ar_ols(x, order, order_max, demean, intercept)
      case ("mle")
         fit = ar_mle(x, order, order_max, demean, mle_delta)
      case default
         fit%status = r_invalid_input
      end select
   end function ar

   pure function spec_ar(x, order, order_max, n_frequency, demean, frequency, method, &
                         intercept, variance_method, mle_delta) result(estimate)
      !! Fits a Yule-Walker AR model and evaluates its spectral density.
      real(dp), intent(in) :: x(:) !! Finite univariate observations.
      integer, intent(in), optional :: order !! Fixed AR order; omission selects an order by AIC.
      integer, intent(in), optional :: order_max !! Largest AIC candidate order.
      integer, intent(in), optional :: n_frequency !! Number of frequencies from zero to Nyquist.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to true.
      real(dp), intent(in), optional :: frequency !! Sampling frequency; defaults to one.
      character(len=*), intent(in), optional :: method !! `yule-walker`, `burg`, `ols`, or `mle`.
      logical, intent(in), optional :: intercept !! OLS intercept; defaults to `demean`.
      integer, intent(in), optional :: variance_method !! Burg variance method one or two.
      real(dp), intent(in), optional :: mle_delta !! MLE steady-state tolerance; defaults to zero.
      type(spectrum_result_t) :: estimate
      type(ar_fit_t) :: fit
      real(dp) :: angle, cosine_sum, sampling_frequency, sine_sum
      integer :: i, j, frequency_count

      sampling_frequency = optval(frequency, 1.0_dp)
      frequency_count = optval(n_frequency, 500)
      if (sampling_frequency <= 0.0_dp .or. frequency_count < 2) then
         estimate%status = r_invalid_input
         return
      end if
      fit = ar(x, order, order_max, demean, method, intercept, variance_method, mle_delta)
      if (fit%status /= r_ok) then
         estimate%status = fit%status
         return
      end if
      allocate (estimate%frequency(frequency_count), estimate%spectrum(frequency_count))
      allocate (estimate%kernel(1), source=1.0_dp)
      do i = 1, frequency_count
         estimate%frequency(i) = 0.5_dp*sampling_frequency*real(i - 1, dp)/ &
                                 real(frequency_count - 1, dp)
         cosine_sum = 0.0_dp
         sine_sum = 0.0_dp
         do j = 1, fit%order
            angle = 2.0_dp*acos(-1.0_dp)*estimate%frequency(i)*real(j, dp)/ &
                    sampling_frequency
            cosine_sum = cosine_sum + fit%coefficients(j)*cos(angle)
            sine_sum = sine_sum + fit%coefficients(j)*sin(angle)
         end do
         estimate%spectrum(i) = fit%variance/ &
            (sampling_frequency*((1.0_dp - cosine_sum)**2 + sine_sum**2))
      end do
      estimate%n_used = size(x)
      estimate%n_original = size(x)
      estimate%order = fit%order
      estimate%demeaned = fit%demeaned
   end function spec_ar

   pure subroutine fit_mle_order(x, order, estimate_mean, delta, coefficients, mean_value, &
                                 variance, log_likelihood, status)
      !! Optimizes one stationary autoregressive order in reflection-coefficient coordinates.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: order !! Fixed autoregressive order.
      logical, intent(in) :: estimate_mean !! Whether to estimate the stationary mean.
      real(dp), intent(in) :: delta !! Steady-state likelihood switch tolerance.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Maximizing AR coefficients.
      real(dp), intent(out) :: mean_value !! Maximizing stationary mean.
      real(dp), intent(out) :: variance !! Maximizing innovation variance.
      real(dp), intent(out) :: log_likelihood !! Maximized Gaussian log likelihood.
      integer, intent(out) :: status !! Completion status from `r_status`.
      real(dp), allocatable :: ar_table(:, :), autocovariance(:), partial(:)
      real(dp), allocatable :: variances(:)
      integer :: local_status

      status = r_ok
      if (order == 0) then
         allocate (coefficients(0))
         call stationary_ar_likelihood(x, coefficients, estimate_mean, delta, &
                                       mean_value, variance, log_likelihood, status)
         return
      end if
      call r_autocovariance(x, autocovariance, lag_max=order, &
                            demean=estimate_mean, status=local_status)
      if (local_status /= r_ok) then
         status = local_status
         return
      end if
      call r_durbin_levinson(autocovariance, partial, local_status, ar_table, variances)
      if (local_status /= r_ok) then
         status = local_status
         return
      end if
      partial = max(-0.98_dp, min(0.98_dp, partial))
      call nelder_mead_mle(x, estimate_mean, delta, partial, local_status)
      if (local_status /= r_ok) then
         status = local_status
         return
      end if
      call reflection_to_ar(partial, coefficients)
      call stationary_ar_likelihood(x, coefficients, estimate_mean, delta, mean_value, variance, &
                                    log_likelihood, status)
   end subroutine fit_mle_order

   pure subroutine nelder_mead_mle(x, estimate_mean, delta, parameters, status)
      !! Minimizes the concentrated likelihood with a bounded Nelder-Mead simplex.
      real(dp), intent(in) :: x(:) !! Finite observations.
      logical, intent(in) :: estimate_mean !! Whether the stationary mean is estimated.
      real(dp), intent(in) :: delta !! Steady-state likelihood switch tolerance.
      real(dp), intent(inout) :: parameters(:) !! Initial and optimized reflection coefficients.
      integer, intent(out) :: status !! Completion status from `r_status`.
      real(dp), allocatable :: centroid(:), contracted(:), expanded(:), reflected(:)
      real(dp), allocatable :: simplex(:, :), values(:)
      real(dp) :: contracted_value, expanded_value, reflected_value
      integer :: i, iteration, local_status, n

      n = size(parameters)
      allocate (simplex(n, n + 1), values(n + 1))
      simplex(:, 1) = parameters
      do i = 1, n
         simplex(:, i + 1) = parameters
         simplex(i, i + 1) = max(-0.999_dp, min(0.999_dp, parameters(i) + 0.1_dp))
         if (simplex(i, i + 1) == parameters(i)) then
            simplex(i, i + 1) = parameters(i) - 0.1_dp
         end if
      end do
      do i = 1, n + 1
         call mle_objective(x, simplex(:, i), estimate_mean, delta, values(i), local_status)
         if (local_status /= r_ok) then
            status = local_status
            return
         end if
      end do
      allocate (centroid(n), contracted(n), expanded(n), reflected(n))
      do iteration = 1, 2000
         call sort_simplex(simplex, values)
         if (maxval(abs(simplex(:, 2:) - spread(simplex(:, 1), 2, n))) < 1.0e-9_dp .and. &
             maxval(abs(values(2:) - values(1))) < 1.0e-10_dp) exit
         centroid = sum(simplex(:, 1:n), dim=2)/real(n, dp)
         reflected = bounded_reflections(centroid + (centroid - simplex(:, n + 1)))
         call mle_objective(x, reflected, estimate_mean, delta, reflected_value, local_status)
         if (reflected_value < values(1)) then
            expanded = bounded_reflections(centroid + 2.0_dp*(reflected - centroid))
            call mle_objective(x, expanded, estimate_mean, delta, expanded_value, local_status)
            if (expanded_value < reflected_value) then
               simplex(:, n + 1) = expanded
               values(n + 1) = expanded_value
            else
               simplex(:, n + 1) = reflected
               values(n + 1) = reflected_value
            end if
         else if (reflected_value < values(n)) then
            simplex(:, n + 1) = reflected
            values(n + 1) = reflected_value
         else
            if (reflected_value < values(n + 1)) then
               contracted = bounded_reflections(centroid + 0.5_dp*(reflected - centroid))
            else
               contracted = bounded_reflections(centroid + 0.5_dp*(simplex(:, n + 1) - centroid))
            end if
            call mle_objective(x, contracted, estimate_mean, delta, contracted_value, local_status)
            if (contracted_value < min(reflected_value, values(n + 1))) then
               simplex(:, n + 1) = contracted
               values(n + 1) = contracted_value
            else
               do i = 2, n + 1
                  simplex(:, i) = bounded_reflections(simplex(:, 1) + &
                                                       0.5_dp*(simplex(:, i) - simplex(:, 1)))
                  call mle_objective(x, simplex(:, i), estimate_mean, delta, &
                                     values(i), local_status)
               end do
            end if
         end if
      end do
      call sort_simplex(simplex, values)
      parameters = simplex(:, 1)
      status = r_ok
   end subroutine nelder_mead_mle

   pure subroutine sort_simplex(simplex, values)
      !! Sorts simplex vertices in ascending objective order.
      real(dp), intent(inout) :: simplex(:, :) !! Parameter vertices arranged by column.
      real(dp), intent(inout) :: values(:) !! Objective values corresponding to vertices.
      real(dp), allocatable :: temporary(:)
      real(dp) :: temporary_value
      integer :: i, j

      allocate (temporary(size(simplex, 1)))
      do i = 2, size(values)
         j = i
         do while (j > 1)
            if (values(j) >= values(j - 1)) exit
            temporary = simplex(:, j)
            simplex(:, j) = simplex(:, j - 1)
            simplex(:, j - 1) = temporary
            temporary_value = values(j)
            values(j) = values(j - 1)
            values(j - 1) = temporary_value
            j = j - 1
         end do
      end do
   end subroutine sort_simplex

   pure function bounded_reflections(parameters) result(bounded)
      !! Restricts reflection coefficients to the open stationarity interval.
      real(dp), intent(in) :: parameters(:) !! Candidate reflection coefficients.
      real(dp) :: bounded(size(parameters))

      bounded = max(-0.999_dp, min(0.999_dp, parameters))
   end function bounded_reflections

   pure subroutine mle_objective(x, partial, estimate_mean, delta, objective, status)
      !! Evaluates the concentrated stationary Gaussian objective in reflection coordinates.
      real(dp), intent(in) :: x(:) !! Finite observations.
      real(dp), intent(in) :: partial(:) !! Reflection coefficients within the unit interval.
      logical, intent(in) :: estimate_mean !! Whether to estimate the stationary mean.
      real(dp), intent(in) :: delta !! Steady-state likelihood switch tolerance.
      real(dp), intent(out) :: objective !! Negative maximized log likelihood.
      integer, intent(out) :: status !! Completion status from `r_status`.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: log_likelihood, mean_value, variance

      call reflection_to_ar(partial, coefficients)
      call stationary_ar_likelihood(x, coefficients, estimate_mean, delta, mean_value, variance, &
                                    log_likelihood, status)
      objective = -log_likelihood
   end subroutine mle_objective

   pure subroutine stationary_ar_likelihood(x, coefficients, estimate_mean, delta, mean_value, &
                                            variance, log_likelihood, status)
      !! Evaluates the exact stationary Gaussian likelihood for fixed AR coefficients.
      real(dp), intent(in) :: x(:) !! Finite observations.
      real(dp), intent(in) :: coefficients(:) !! Stationary AR coefficients.
      logical, intent(in) :: estimate_mean !! Whether to estimate the stationary mean.
      real(dp), intent(in) :: delta !! Switch to steady-state recursion below this variance excess.
      real(dp), intent(out) :: mean_value !! Generalized least-squares mean.
      real(dp), intent(out) :: variance !! Concentrated innovation variance.
      real(dp), intent(out) :: log_likelihood !! Maximized Gaussian log likelihood.
      integer, intent(out) :: status !! Completion status from `r_status`.
      real(dp), allocatable :: autocorrelation(:), covariance(:, :), innovation_one(:)
      real(dp), allocatable :: innovation_x(:), moving_average(:), state_one(:), state_x(:)
      real(dp), allocatable :: transition(:, :)
      real(dp) :: denominator, forecast_variance, gamma_zero, sum_log_variance
      real(dp) :: weighted_cross, weighted_one
      integer :: i, j, order
      logical :: steady_state

      allocate (moving_average(0))
      call r_arma_autocorrelation(coefficients, moving_average, size(x) - 1, &
                                  autocorrelation, status)
      if (status /= r_ok) return
      denominator = 1.0_dp
      if (size(coefficients) > 0) then
         denominator = denominator - dot_product(coefficients, &
            autocorrelation(1:size(coefficients)))
      end if
      if (denominator <= tiny(1.0_dp)) then
         status = r_invalid_input
         return
      end if
      gamma_zero = 1.0_dp/denominator
      order = size(coefficients)
      allocate (innovation_x(size(x)), innovation_one(size(x)))
      if (order == 0) then
         innovation_x = x
         innovation_one = 1.0_dp
         sum_log_variance = 0.0_dp
      else
         allocate (covariance(order, order), transition(order, order), source=0.0_dp)
         allocate (state_x(order), state_one(order), source=0.0_dp)
         do j = 1, order
            do i = 1, order
               covariance(i, j) = gamma_zero*autocorrelation(abs(i - j))
            end do
         end do
         transition(1, :) = coefficients
         do i = 2, order
            transition(i, i - 1) = 1.0_dp
         end do
         sum_log_variance = 0.0_dp
         steady_state = .false.
         do i = 1, size(x)
            forecast_variance = covariance(1, 1)
            if (steady_state .or. (delta > 0.0_dp .and. forecast_variance - 1.0_dp < delta)) then
               steady_state = .true.
               innovation_x(i) = x(i)
               innovation_one(i) = 1.0_dp
               do j = 1, min(order, i - 1)
                  innovation_x(i) = innovation_x(i) - coefficients(j)*x(i - j)
                  innovation_one(i) = innovation_one(i) - coefficients(j)
               end do
            else
               innovation_x(i) = (x(i) - state_x(1))/sqrt(forecast_variance)
               innovation_one(i) = (1.0_dp - state_one(1))/sqrt(forecast_variance)
               state_x = state_x + covariance(:, 1)*(x(i) - state_x(1))/forecast_variance
               state_one = state_one + covariance(:, 1)*(1.0_dp - state_one(1))/forecast_variance
               covariance = covariance - spread(covariance(:, 1), 2, order)* &
                            spread(covariance(1, :), 1, order)/forecast_variance
               state_x = matmul(transition, state_x)
               state_one = matmul(transition, state_one)
               covariance = matmul(transition, matmul(covariance, transpose(transition)))
               covariance(1, 1) = covariance(1, 1) + 1.0_dp
               sum_log_variance = sum_log_variance + log(forecast_variance)
            end if
         end do
      end if
      mean_value = 0.0_dp
      if (estimate_mean) then
         weighted_cross = dot_product(innovation_x, innovation_one)
         weighted_one = dot_product(innovation_one, innovation_one)
         mean_value = weighted_cross/weighted_one
      end if
      variance = sum((innovation_x - mean_value*innovation_one)**2)/real(size(x), dp)
      if (variance <= tiny(1.0_dp)) then
         status = r_invalid_input
         return
      end if
      log_likelihood = -0.5_dp*(real(size(x), dp)*(log(2.0_dp*acos(-1.0_dp)*variance) + 1.0_dp) + &
                                 sum_log_variance)
      status = r_ok
   end subroutine stationary_ar_likelihood

   pure subroutine reflection_to_ar(partial, coefficients)
      !! Converts reflection coefficients to stationary autoregressive coefficients.
      real(dp), intent(in) :: partial(:) !! Reflection coefficients in increasing order.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Stationary AR coefficients.
      real(dp), allocatable :: previous(:)
      integer :: i, j

      allocate (coefficients(size(partial)), source=0.0_dp)
      do i = 1, size(partial)
         if (i > 1) previous = coefficients
         coefficients(i) = partial(i)
         do j = 1, i - 1
            coefficients(j) = previous(j) - partial(i)*previous(i - j)
         end do
      end do
   end subroutine reflection_to_ar

   pure subroutine calculate_residuals(x, coefficients, mean_value, intercept, residuals)
      !! Calculates fitted autoregressive residuals with leading NaNs.
      real(dp), intent(in) :: x(:) !! Original observations.
      real(dp), intent(in) :: coefficients(:) !! AR coefficients in lag order.
      real(dp), intent(in) :: mean_value !! Mean removed before fitting.
      real(dp), intent(in) :: intercept !! Intercept fitted to the centered observations.
      real(dp), allocatable, intent(out) :: residuals(:) !! Residual sequence with leading NaNs.
      integer :: i, j

      allocate (residuals(size(x)))
      residuals = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = size(coefficients) + 1, size(x)
         residuals(i) = x(i) - mean_value - intercept
         do j = 1, size(coefficients)
            residuals(i) = residuals(i) - coefficients(j)*(x(i - j) - mean_value)
         end do
      end do
   end subroutine calculate_residuals

   pure subroutine calculate_ols_intercept(y, order, include_intercept, coefficients, intercept)
      !! Recomputes the selected OLS intercept without retaining every candidate fit.
      real(dp), intent(in) :: y(:) !! Demeaned or original observations used for fitting.
      integer, intent(in) :: order !! Selected autoregressive order.
      logical, intent(in) :: include_intercept !! Whether an intercept is fitted.
      real(dp), intent(in) :: coefficients(:) !! Selected autoregressive coefficients.
      real(dp), intent(out) :: intercept !! Selected regression intercept.
      integer :: i

      intercept = 0.0_dp
      if (.not. include_intercept) return
      intercept = sum(y(order + 1:))
      do i = 1, order
         intercept = intercept - coefficients(i)*sum(y(order + 1 - i:size(y) - i))
      end do
      intercept = intercept/real(size(y) - order, dp)
   end subroutine calculate_ols_intercept

   pure function lowercase(text) result(lower)
      !! Converts ASCII uppercase letters to lowercase.
      character(len=*), intent(in) :: text !! Input text.
      character(len=len(text)) :: lower
      integer :: code, i

      lower = text
      do i = 1, len(text)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

end module r_stats_autoregression
