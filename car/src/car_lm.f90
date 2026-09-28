module car_lm
   use car_kinds, only: dp
   use car_utils, only: f_survival, least_squares_residuals, matrix_inverse, normal_quantile, &
                        student_t_quantile, student_t_two_sided_p
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   implicit none
   private

   public :: confint_lm
   public :: lm_bootstrap_case
   public :: lm_refit_rows
   public :: lm_summary
   public :: predict_lm

contains

   pure subroutine lm_summary(x, y, coefficients, standard_errors, t_values, p_values, sigma, &
                              r_squared, adjusted_r_squared, f_statistic, f_p_value, residual_df, &
                              status, weights, has_intercept, covariance)
      real(dp), intent(in) :: x(:, :) !! Design matrix with observations in rows and coefficients in columns.
      real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Least-squares coefficient estimates in design-column order.
      real(dp), allocatable, intent(out) :: standard_errors(:) !! Standard errors from the classical least-squares covariance matrix.
      real(dp), allocatable, intent(out) :: t_values(:) !! Coefficient t statistics, estimate divided by standard error.
      real(dp), allocatable, intent(out) :: p_values(:) !! Two-sided Student-t p-values for the coefficient tests.
      real(dp), intent(out) :: sigma !! Residual standard error sqrt(RSS / residual_df) on the weighted scale.
      real(dp), intent(out) :: r_squared !! Coefficient of determination using weighted or unweighted sums of squares.
      real(dp), intent(out) :: adjusted_r_squared !! Degrees-of-freedom-adjusted coefficient of determination.
      real(dp), intent(out) :: f_statistic !! Overall regression F statistic excluding the intercept when present.
      real(dp), intent(out) :: f_p_value !! Upper-tail probability for the overall F statistic.
      integer, intent(out) :: residual_df !! Residual degrees of freedom n - p for the full-rank fit.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, weights, or singular normal equations.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights; omitted means ordinary least squares.
      logical, intent(in), optional :: has_intercept !! Whether the model includes an intercept; inferred from a constant first column when omitted.
      real(dp), allocatable, intent(out), optional :: covariance(:, :) !! Optional classical covariance matrix of the coefficient estimates.
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: xtwx(:, :)
      real(dp), allocatable :: xtwx_inverse(:, :)
      real(dp) :: f_num
      real(dp) :: mean_y
      real(dp) :: rss
      real(dp) :: sum_w
      real(dp) :: tss
      real(dp) :: weight
      integer :: i
      integer :: inv_status
      integer :: n
      integer :: numerator_df
      integer :: p
      logical :: intercept

      n = size(x, 1)
      p = size(x, 2)
      allocate(coefficients(p), standard_errors(p), t_values(p), p_values(p))
      coefficients = 0.0_dp
      standard_errors = 0.0_dp
      t_values = 0.0_dp
      p_values = 1.0_dp
      sigma = 0.0_dp
      r_squared = 0.0_dp
      adjusted_r_squared = 0.0_dp
      f_statistic = 0.0_dp
      f_p_value = 1.0_dp
      residual_df = n - p
      status = 0

      if (n <= 0 .or. p <= 0 .or. size(y) /= n .or. residual_df <= 0) then
         status = 1
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = 2
            return
         end if
         if (sum(weights) <= 0.0_dp) then
            status = 2
            return
         end if
      end if

      allocate(residuals(n), xtwx(p, p), xtwx_inverse(p, p))
      if (present(weights)) then
         call least_squares_residuals(x, y, coefficients, residuals, status, weights)
      else
         call least_squares_residuals(x, y, coefficients, residuals, status)
      end if
      if (status /= 0) return

      xtwx = 0.0_dp
      do i = 1, n
         if (present(weights)) then
            weight = weights(i)
         else
            weight = 1.0_dp
         end if
         xtwx = xtwx + weight * outer_product(x(i, :), x(i, :))
      end do
      call matrix_inverse(xtwx, xtwx_inverse, inv_status)
      if (inv_status /= 0) then
         status = 3
         return
      end if

      rss = dot_product(residuals, residuals)
      sigma = sqrt(rss / real(residual_df, dp))
      do i = 1, p
         standard_errors(i) = sigma * sqrt(max(0.0_dp, xtwx_inverse(i, i)))
         if (standard_errors(i) > 0.0_dp) then
            t_values(i) = coefficients(i) / standard_errors(i)
            p_values(i) = student_t_two_sided_p(t_values(i), residual_df)
         else if (abs(coefficients(i)) > 0.0_dp) then
            t_values(i) = sign(huge(1.0_dp), coefficients(i))
            p_values(i) = 0.0_dp
         end if
      end do
      if (present(covariance)) covariance = sigma * sigma * xtwx_inverse

      if (present(has_intercept)) then
         intercept = has_intercept
      else
         intercept = is_constant_column(x(:, 1))
      end if

      if (present(weights)) then
         sum_w = sum(weights)
         if (intercept) then
            mean_y = dot_product(weights, y) / sum_w
            tss = dot_product(weights, (y - mean_y)**2)
         else
            tss = dot_product(weights, y * y)
         end if
      else
         if (intercept) then
            mean_y = sum(y) / real(n, dp)
            tss = sum((y - mean_y)**2)
         else
            tss = sum(y * y)
         end if
      end if

      if (tss > 0.0_dp) then
         r_squared = 1.0_dp - rss / tss
      else
         r_squared = 0.0_dp
      end if
      if (intercept) then
         adjusted_r_squared = 1.0_dp - (1.0_dp - r_squared) * real(n - 1, dp) / real(residual_df, dp)
         numerator_df = p - 1
      else
         adjusted_r_squared = 1.0_dp - (1.0_dp - r_squared) * real(n, dp) / real(residual_df, dp)
         numerator_df = p
      end if

      if (numerator_df > 0 .and. rss > 0.0_dp) then
         f_num = max(0.0_dp, tss - rss) / real(numerator_df, dp)
         f_statistic = f_num / (rss / real(residual_df, dp))
         f_p_value = f_survival(f_statistic, numerator_df, residual_df)
      end if
   end subroutine lm_summary

   pure subroutine confint_lm(coefficients, covariance, residual_df, level, lower, upper, status)
      real(dp), intent(in) :: coefficients(:) !! Estimated regression coefficients.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix for coefficients in the same order.
      integer, intent(in) :: residual_df !! Positive residual degrees of freedom for Student-t intervals.
      real(dp), intent(in) :: level !! Central confidence level strictly between zero and one.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower confidence limit for each coefficient.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper confidence limit for each coefficient.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, level, or degrees of freedom.
      real(dp) :: critical
      real(dp) :: se
      integer :: i
      integer :: p
      integer :: quantile_status

      p = size(coefficients)
      allocate(lower(p), upper(p))
      lower = 0.0_dp
      upper = 0.0_dp
      status = 0
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p .or. p <= 0) then
         status = 1
         return
      end if
      if (residual_df <= 0 .or. level <= 0.0_dp .or. level >= 1.0_dp) then
         status = 2
         return
      end if
      if (any([(covariance(i, i) < 0.0_dp, i = 1, p)])) then
         status = 3
         return
      end if

      call student_t_quantile(0.5_dp * (1.0_dp + level), residual_df, critical, quantile_status)
      if (quantile_status /= 0) then
         status = 4
         return
      end if
      do i = 1, p
         se = sqrt(max(0.0_dp, covariance(i, i)))
         lower(i) = coefficients(i) - critical * se
         upper(i) = coefficients(i) + critical * se
      end do
   end subroutine confint_lm

   pure subroutine predict_lm(coefficients, covariance, x_new, fit, standard_error, lower, upper, status, &
                              level, residual_df, residual_variance, prediction_interval, offset)
      real(dp), intent(in) :: coefficients(:) !! Regression coefficients in the same order as columns of x_new.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix of coefficients.
      real(dp), intent(in) :: x_new(:, :) !! New design matrix with one prediction case per row.
      real(dp), allocatable, intent(out) :: fit(:) !! Fitted response for each row of x_new, including any supplied offset.
      real(dp), allocatable, intent(out) :: standard_error(:) !! Standard error of each fitted mean response.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower confidence or prediction interval limit.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper confidence or prediction interval limit.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid shapes, options, or variances.
      real(dp), intent(in), optional :: level !! Central interval level; default is 0.95.
      integer, intent(in), optional :: residual_df !! Positive residual degrees of freedom; omitted uses a normal critical value.
      real(dp), intent(in), optional :: residual_variance !! Nonnegative future-response variance added for prediction intervals.
      logical, intent(in), optional :: prediction_interval !! True requests a prediction interval; default is confidence interval.
      real(dp), intent(in), optional :: offset(:) !! Optional additive offset with one value per prediction row.
      real(dp) :: alpha
      real(dp) :: critical
      real(dp) :: interval_variance
      real(dp) :: level_value
      real(dp) :: mean_variance
      integer :: i
      integer :: n
      integer :: p
      integer :: quantile_status
      logical :: prediction

      n = size(x_new, 1)
      p = size(x_new, 2)
      allocate(fit(n), standard_error(n), lower(n), upper(n))
      fit = 0.0_dp
      standard_error = 0.0_dp
      lower = 0.0_dp
      upper = 0.0_dp
      status = 0
      if (size(coefficients) /= p .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p) then
         status = 1
         return
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            status = 2
            return
         end if
      end if

      level_value = 0.95_dp
      if (present(level)) level_value = level
      if (level_value <= 0.0_dp .or. level_value >= 1.0_dp) then
         status = 3
         return
      end if
      alpha = 0.5_dp * (1.0_dp + level_value)
      if (present(residual_df)) then
         if (residual_df <= 0) then
            status = 4
            return
         end if
         call student_t_quantile(alpha, residual_df, critical, quantile_status)
         if (quantile_status /= 0) then
            status = 4
            return
         end if
      else
         critical = normal_quantile(alpha)
      end if

      prediction = .false.
      if (present(prediction_interval)) prediction = prediction_interval
      if (prediction) then
         if (.not. present(residual_variance)) then
            status = 5
            return
         end if
         if (residual_variance < 0.0_dp) then
            status = 5
            return
         end if
      end if

      fit = matmul(x_new, coefficients)
      if (present(offset)) fit = fit + offset
      do i = 1, n
         mean_variance = dot_product(x_new(i, :), matmul(covariance, x_new(i, :)))
         if (mean_variance < -100.0_dp * epsilon(1.0_dp)) then
            status = 6
            return
         end if
         mean_variance = max(0.0_dp, mean_variance)
         standard_error(i) = sqrt(mean_variance)
         interval_variance = mean_variance
         if (prediction) interval_variance = interval_variance + residual_variance
         lower(i) = fit(i) - critical * sqrt(interval_variance)
         upper(i) = fit(i) + critical * sqrt(interval_variance)
      end do
   end subroutine predict_lm

   pure subroutine lm_refit_rows(x, y, indices, coefficients, status, weights)
      real(dp), intent(in) :: x(:, :) !! Original design matrix with observations in rows and coefficients in columns.
      real(dp), intent(in) :: y(:) !! Original response vector with one value per row of x.
      integer, intent(in) :: indices(:) !! One-based row indices defining the resampled or subsetted refit, with repetition allowed.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Least-squares coefficients fitted to the selected rows.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid indices, weights, or a singular refit.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative original case weights selected in parallel with rows.
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: selected_weights(:)
      real(dp), allocatable :: x_selected(:, :)
      real(dp), allocatable :: y_selected(:)
      integer :: i
      integer :: n
      integer :: n_selected
      integer :: p

      n = size(x, 1)
      p = size(x, 2)
      n_selected = size(indices)
      allocate(coefficients(p))
      coefficients = ieee_value(0.0_dp, ieee_quiet_nan)
      status = 0
      if (n <= 0 .or. p <= 0 .or. size(y) /= n .or. n_selected <= 0) then
         status = 1
         return
      end if
      if (any(indices < 1) .or. any(indices > n)) then
         status = 2
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = 3
            return
         end if
      end if

      allocate(x_selected(n_selected, p), y_selected(n_selected), residuals(n_selected))
      if (present(weights)) allocate(selected_weights(n_selected))
      do i = 1, n_selected
         x_selected(i, :) = x(indices(i), :)
         y_selected(i) = y(indices(i))
         if (present(weights)) selected_weights(i) = weights(indices(i))
      end do
      if (present(weights)) then
         call least_squares_residuals(x_selected, y_selected, coefficients, residuals, status, selected_weights)
      else
         call least_squares_residuals(x_selected, y_selected, coefficients, residuals, status)
      end if
   end subroutine lm_refit_rows

   pure subroutine lm_bootstrap_case(x, y, sample_indices, estimates, status, weights)
      real(dp), intent(in) :: x(:, :) !! Original design matrix with observations in rows and coefficients in columns.
      real(dp), intent(in) :: y(:) !! Original response vector with one value per row of x.
      integer, intent(in) :: sample_indices(:, :) !! Bootstrap row indices; each column is one size-n resample with replacement.
      real(dp), allocatable, intent(out) :: estimates(:, :) !! Bootstrap coefficient estimates with replicates in rows and coefficients in columns.
      integer, intent(out) :: status !! Zero on success; negative for invalid input, or the one-based replicate number of a singular fit.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights resampled with x and y.
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: xb(:, :)
      real(dp), allocatable :: wb(:)
      real(dp), allocatable :: yb(:)
      integer :: b
      integer :: fit_status
      integer :: i
      integer :: n
      integer :: p
      integer :: replicates

      n = size(x, 1)
      p = size(x, 2)
      replicates = size(sample_indices, 2)
      allocate(estimates(replicates, p))
      estimates = ieee_value(0.0_dp, ieee_quiet_nan)
      status = 0
      if (n <= 0 .or. p <= 0 .or. size(y) /= n .or. size(sample_indices, 1) /= n .or. replicates <= 0) then
         status = -1
         return
      end if
      if (any(sample_indices < 1) .or. any(sample_indices > n)) then
         status = -2
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = -3
            return
         end if
      end if

      allocate(beta(p), residuals(n), xb(n, p), yb(n))
      if (present(weights)) allocate(wb(n))
      do b = 1, replicates
         do i = 1, n
            xb(i, :) = x(sample_indices(i, b), :)
            yb(i) = y(sample_indices(i, b))
            if (present(weights)) wb(i) = weights(sample_indices(i, b))
         end do
         if (present(weights)) then
            call least_squares_residuals(xb, yb, beta, residuals, fit_status, wb)
         else
            call least_squares_residuals(xb, yb, beta, residuals, fit_status)
         end if
         if (fit_status /= 0) then
            status = b
            return
         end if
         estimates(b, :) = beta
      end do
   end subroutine lm_bootstrap_case

   pure function outer_product(a, b) result(product)
      real(dp), intent(in) :: a(:) !! Left vector in the outer product.
      real(dp), intent(in) :: b(:) !! Right vector in the outer product.
      real(dp) :: product(size(a), size(b))
      integer :: i

      do i = 1, size(a)
         product(i, :) = a(i) * b
      end do
   end function outer_product

   pure function is_constant_column(x) result(is_constant)
      real(dp), intent(in) :: x(:) !! Candidate design column tested for numerical constancy.
      logical :: is_constant

      if (size(x) == 0) then
         is_constant = .false.
      else
         is_constant = maxval(abs(x - x(1))) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(1)))
      end if
   end function is_constant_column

end module car_lm
