module car_diagnostics
   use car_kinds, only: dp
   use car_utils, only: chi_square_survival, f_survival, least_squares_residuals, mean_valid, median_valid, &
                        student_t_two_sided_p
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private

   public :: durbin_watson
   public :: levene_test
   public :: ncv_test
   public :: outlier_bonferroni

contains

   subroutine durbin_watson(residuals, max_lag, dw, autocorrelation, status)
      real(dp), intent(in) :: residuals(:) !! Ordered model residuals with no missing values.
      integer, intent(in) :: max_lag !! Largest positive lag, not exceeding size(residuals)-1.
      real(dp), allocatable, intent(out) :: dw(:) !! Generalized Durbin-Watson statistics for lags 1:max_lag.
      real(dp), allocatable, intent(out), optional :: autocorrelation(:) !! Residual cross-products divided by total residual sum of squares.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid lag, NaN residual, or zero denominator.
      real(dp) :: denominator
      integer :: lag, n

      status = 0
      n = size(residuals)
      if (max_lag < 1 .or. max_lag >= n .or. any(ieee_is_nan(residuals))) then
         status = 1
         allocate(dw(0))
         if (present(autocorrelation)) allocate(autocorrelation(0))
         return
      end if
      denominator = sum(residuals**2)
      if (denominator <= tiny(1.0_dp)) then
         status = 2
         allocate(dw(0))
         if (present(autocorrelation)) allocate(autocorrelation(0))
         return
      end if
      allocate(dw(max_lag))
      if (present(autocorrelation)) allocate(autocorrelation(max_lag))
      do lag = 1, max_lag
         dw(lag) = sum((residuals(lag + 1:n) - residuals(1:n - lag))**2) / denominator
         if (present(autocorrelation)) then
            autocorrelation(lag) = sum(residuals(lag + 1:n) * residuals(1:n - lag)) / denominator
         end if
      end do
   end subroutine durbin_watson

   subroutine levene_test(y, group, statistic, p_value, df1, df2, status, center)
      real(dp), intent(in) :: y(:) !! Numeric response observations.
      integer, intent(in) :: group(:) !! Positive integer group labels, one per observation.
      real(dp), intent(out) :: statistic !! Levene one-way ANOVA F statistic for absolute center deviations.
      real(dp), intent(out) :: p_value !! Upper-tail F probability for the Levene statistic.
      integer, intent(out) :: df1 !! Numerator degrees of freedom, number of nonempty groups minus one.
      integer, intent(out) :: df2 !! Denominator degrees of freedom, valid observations minus number of groups.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, groups, or degenerate variation.
      character(len=*), intent(in), optional :: center !! Group center: "median" (default) or "mean".
      real(dp), allocatable :: centers(:), deviation(:)
      real(dp) :: grand, ss_between, ss_within, center_value, group_mean
      integer :: g, k, n, ng, count_g, n_valid, center_status, used
      character(len=8) :: center_kind
      logical, allocatable :: valid(:), in_group(:)

      statistic = 0.0_dp
      p_value = 1.0_dp
      df1 = 0
      df2 = 0
      status = 0
      n = size(y)
      if (size(group) /= n .or. n < 2 .or. any(group <= 0)) then
         status = 1
         return
      end if
      center_kind = "median"
      if (present(center)) center_kind = lowercase(adjustl(center))
      if (trim(center_kind) /= "median" .and. trim(center_kind) /= "mean") then
         status = 2
         return
      end if

      valid = .not. ieee_is_nan(y)
      n_valid = count(valid)
      if (n_valid < 2) then
         status = 3
         return
      end if
      ng = maxval(group, mask=valid)
      allocate(centers(ng), deviation(n), in_group(n))
      centers = 0.0_dp
      deviation = 0.0_dp
      used = 0
      do g = 1, ng
         in_group = valid .and. group == g
         count_g = count(in_group)
         if (count_g == 0) cycle
         used = used + 1
         if (trim(center_kind) == "median") then
            call median_valid(pack(y, in_group), center_value, k, center_status)
         else
            call mean_valid(pack(y, in_group), center_value, k, center_status)
         end if
         if (center_status /= 0) then
            status = 3
            return
         end if
         centers(g) = center_value
      end do
      if (used < 2 .or. n_valid <= used) then
         status = 4
         return
      end if
      do k = 1, n
         if (valid(k)) deviation(k) = abs(y(k) - centers(group(k)))
      end do
      grand = sum(pack(deviation, valid)) / real(n_valid, dp)
      ss_between = 0.0_dp
      ss_within = 0.0_dp
      do g = 1, ng
         in_group = valid .and. group == g
         count_g = count(in_group)
         if (count_g == 0) cycle
         group_mean = sum(pack(deviation, in_group)) / real(count_g, dp)
         ss_between = ss_between + real(count_g, dp) * (group_mean - grand)**2
         ss_within = ss_within + sum((pack(deviation, in_group) - group_mean)**2)
      end do
      df1 = used - 1
      df2 = n_valid - used
      if (ss_within <= tiny(1.0_dp)) then
         status = 5
         return
      end if
      statistic = (ss_between / real(df1, dp)) / (ss_within / real(df2, dp))
      p_value = f_survival(statistic, df1, df2)
   end subroutine levene_test

   subroutine ncv_test(residuals, fitted, residual_df, chi_square, p_value, status)
      real(dp), intent(in) :: residuals(:) !! Pearson or ordinary residuals from the fitted linear model.
      real(dp), intent(in) :: fitted(:) !! Fitted values used as the default variance predictor.
      integer, intent(in) :: residual_df !! Positive residual degrees of freedom from the original model.
      real(dp), intent(out) :: chi_square !! Score-test statistic for nonconstant variance.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability with one degree of freedom.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid data, scale, or regression fit.
      real(dp), allocatable :: x(:, :), u(:), beta(:), fit_residuals(:)
      real(dp) :: s_sq, mean_u, total_ss, error_ss, regression_ss
      integer :: n, fit_status

      chi_square = 0.0_dp
      p_value = 1.0_dp
      status = 0
      n = size(residuals)
      if (size(fitted) /= n .or. n < 3 .or. residual_df <= 0 .or. &
          any(ieee_is_nan(residuals)) .or. any(ieee_is_nan(fitted))) then
         status = 1
         return
      end if
      s_sq = real(residual_df, dp) * sum(residuals**2) / real(n * residual_df, dp)
      if (s_sq <= tiny(1.0_dp)) then
         status = 2
         return
      end if
      allocate(u(n), x(n, 2), beta(2), fit_residuals(n))
      u = residuals**2 / s_sq
      x(:, 1) = 1.0_dp
      x(:, 2) = fitted
      call least_squares_residuals(x, u, beta, fit_residuals, fit_status)
      if (fit_status /= 0) then
         status = 3
         return
      end if
      mean_u = sum(u) / real(n, dp)
      total_ss = sum((u - mean_u)**2)
      error_ss = sum(fit_residuals**2)
      regression_ss = max(0.0_dp, total_ss - error_ss)
      chi_square = regression_ss / 2.0_dp
      p_value = chi_square_survival(chi_square, 1)
   end subroutine ncv_test

   subroutine outlier_bonferroni(rstudent, residual_df, p_value, bonferroni_p, status)
      real(dp), intent(in) :: rstudent(:) !! Externally studentized residuals from a linear-model fit.
      integer, intent(in) :: residual_df !! Degrees of freedom for the externally studentized t reference distribution.
      real(dp), allocatable, intent(out) :: p_value(:) !! Two-sided unadjusted t probabilities for each residual.
      real(dp), allocatable, intent(out) :: bonferroni_p(:) !! Bonferroni-adjusted probabilities, capped at one.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid degrees of freedom or NaN residuals.
      integer :: i, n

      status = 0
      n = size(rstudent)
      if (residual_df <= 0 .or. n == 0 .or. any(ieee_is_nan(rstudent))) then
         status = 1
         allocate(p_value(0), bonferroni_p(0))
         return
      end if
      allocate(p_value(n), bonferroni_p(n))
      do i = 1, n
         p_value(i) = student_t_two_sided_p(rstudent(i), residual_df)
         bonferroni_p(i) = min(1.0_dp, real(n, dp) * p_value(i))
      end do
   end subroutine outlier_bonferroni

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text !! Option text converted to lowercase.
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

end module car_diagnostics
