module car_power_fit
   use car_kinds, only: dp
   use car_transforms, only: bc_power, bcn_power, yj_power
   use car_utils, only: chi_square_quantile, chi_square_survival, least_squares_residuals, sample_variance
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   implicit none
   private

   public :: box_cox_profile
   public :: inv_transform_estimate
   public :: power_transform_fit
   public :: power_transform_test

contains

   subroutine power_transform_fit(x, y, family, lambda, objective, status, weights, gamma, lower, upper)
      real(dp), intent(in) :: x(:, :) !! Design matrix used to residualize the transformed response.
      real(dp), intent(in) :: y(:) !! Univariate response transformed toward normality.
      character(len=*), intent(in) :: family !! Transformation family: bcPower, yjPower, or bcnPower.
      real(dp), intent(out) :: lambda !! Estimated power parameter minimizing the car likelihood objective.
      real(dp), intent(out) :: objective !! Objective value at lambda; smaller values are preferred.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid input, transform, fit, or optimization failure.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights.
      real(dp), intent(in), optional :: gamma !! BCN gamma parameter; required only for bcnPower and defaulted to zero.
      real(dp), intent(in), optional :: lower !! Lower lambda search bound; default -3.
      real(dp), intent(in), optional :: upper !! Upper lambda search bound; default +3.
      real(dp) :: lo, hi
      integer :: opt_status

      lambda = ieee_value(0.0_dp, ieee_quiet_nan)
      objective = lambda
      status = validate_fit_input(x, y, weights)
      if (status /= 0) return
      lo = -3.0_dp
      hi = 3.0_dp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      call minimize_transform_fit(x, y, family, .true., weights, gamma, lo, hi, 1.0e-7_dp, &
                                  lambda, objective, opt_status)
      if (opt_status /= 0) then
         status = 10 + opt_status
      else
         status = 0
      end if

   end subroutine power_transform_fit

   subroutine power_transform_test(x, y, family, fitted_lambda, null_lambda, lr_statistic, p_value, status, weights, gamma)
      real(dp), intent(in) :: x(:, :) !! Design matrix used for the fitted power transformation.
      real(dp), intent(in) :: y(:) !! Univariate response used for the fitted power transformation.
      character(len=*), intent(in) :: family !! Transformation family: bcPower, yjPower, or bcnPower.
      real(dp), intent(in) :: fitted_lambda !! Estimated lambda defining the unconstrained objective value.
      real(dp), intent(in) :: null_lambda !! Lambda value to test by a one-degree-of-freedom likelihood-ratio test.
      real(dp), intent(out) :: lr_statistic !! Twice the increase in objective at the null lambda.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability with one degree of freedom.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid input or transformation/fit failure.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights.
      real(dp), intent(in), optional :: gamma !! BCN gamma parameter; default zero when the BCN family is selected.
      real(dp) :: fitted_value, null_value

      lr_statistic = 0.0_dp
      p_value = 1.0_dp
      status = validate_fit_input(x, y, weights)
      if (status /= 0) return
      fitted_value = transform_fit_objective(x, y, family, fitted_lambda, .true., weights, gamma)
      null_value = transform_fit_objective(x, y, family, null_lambda, .true., weights, gamma)
      if (.not. ieee_is_finite(fitted_value) .or. .not. ieee_is_finite(null_value)) then
         status = 2
         return
      end if
      lr_statistic = max(0.0_dp, 2.0_dp * (null_value - fitted_value))
      p_value = chi_square_survival(lr_statistic, 1)
   end subroutine power_transform_test

   subroutine box_cox_profile(x, y, lambda_grid, log_likelihood, status, family, weights, gamma)
      real(dp), intent(in) :: x(:, :) !! Design matrix for the regression whose response is transformed.
      real(dp), intent(in) :: y(:) !! Univariate response to transform.
      real(dp), intent(in) :: lambda_grid(:) !! Lambda values at which the profile log-likelihood is evaluated.
      real(dp), allocatable, intent(out) :: log_likelihood(:) !! Profile log-likelihood values, defined up to an additive constant.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or any failed grid evaluation.
      character(len=*), intent(in), optional :: family !! Transformation family; default bcPower.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights.
      real(dp), intent(in), optional :: gamma !! BCN gamma parameter; default zero when selected.
      character(len=16) :: family_name
      real(dp), allocatable :: transformed(:), beta(:), residuals(:)
      real(dp) :: rss
      integer :: i, fit_status, transform_status, n

      status = validate_fit_input(x, y, weights)
      if (status /= 0) then
         allocate(log_likelihood(0))
         return
      end if
      family_name = "bcPower"
      if (present(family)) family_name = trim(family)
      n = size(y)
      allocate(log_likelihood(size(lambda_grid)), transformed(n), beta(size(x, 2)), residuals(n))
      do i = 1, size(lambda_grid)
         call apply_transform(y, family_name, lambda_grid(i), .true., transformed, transform_status, gamma)
         if (transform_status /= 0) then
            status = 2
            log_likelihood = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (present(weights)) then
            call least_squares_residuals(x, transformed, beta, residuals, fit_status, weights)
         else
            call least_squares_residuals(x, transformed, beta, residuals, fit_status)
         end if
         if (fit_status /= 0) then
            status = 3
            log_likelihood = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         rss = sum(residuals**2)
         if (rss <= tiny(1.0_dp)) then
            status = 4
            log_likelihood = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         log_likelihood(i) = -0.5_dp * real(n, dp) * log(rss)
      end do
      status = 0
   end subroutine box_cox_profile

   subroutine inv_transform_estimate(x, y, family, lambda, status, confidence, lower_ci, upper_ci, gamma)
      real(dp), intent(in) :: x(:) !! Predictor values to transform before fitting a straight-line response model.
      real(dp), intent(in) :: y(:) !! Response values paired with x.
      character(len=*), intent(in) :: family !! Transformation family: bcPower, yjPower, or bcnPower.
      real(dp), intent(out) :: lambda !! Lambda minimizing residual sum of squares over [-10, 10].
      integer, intent(out) :: status !! Zero on success; nonzero for invalid data, transformation, fitting, or CI failure.
      real(dp), intent(in), optional :: confidence !! Optional confidence level strictly between zero and one.
      real(dp), intent(out), optional :: lower_ci !! Lower likelihood confidence endpoint when confidence is present.
      real(dp), intent(out), optional :: upper_ci !! Upper likelihood confidence endpoint when confidence is present.
      real(dp), intent(in), optional :: gamma !! BCN gamma parameter; default zero when the BCN family is selected.
      real(dp) :: best_rss, cutoff, quantile, dev0, dummy
      integer :: opt_status, quantile_status, n

      lambda = ieee_value(0.0_dp, ieee_quiet_nan)
      status = 0
      if (size(x) /= size(y) .or. size(x) < 3 .or. any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) then
         status = 1
         return
      end if
      call minimize_inverse_rss(x, y, family, gamma, -10.0_dp, 10.0_dp, 1.0e-7_dp, &
                                lambda, best_rss, opt_status)
      if (opt_status /= 0) then
         status = 2
         return
      end if
      if (.not. present(confidence)) return
      if (.not. present(lower_ci) .or. .not. present(upper_ci)) then
         status = 3
         return
      end if
      if (confidence <= 0.0_dp .or. confidence >= 1.0_dp) then
         status = 4
         return
      end if
      n = size(x)
      call chi_square_quantile(confidence, 1, quantile, quantile_status)
      if (quantile_status /= 0) then
         status = 5
         return
      end if
      dev0 = -real(n, dp) * log(best_rss)
      cutoff = quantile / 2.0_dp
      call minimize_inverse_ci(x, y, family, gamma, dev0, real(n, dp), cutoff, -10.0_dp, lambda, &
                               1.0e-7_dp, lower_ci, dummy, opt_status)
      if (opt_status /= 0) then
         status = 6
         return
      end if
      call minimize_inverse_ci(x, y, family, gamma, dev0, real(n, dp), cutoff, lambda, 10.0_dp, &
                               1.0e-7_dp, upper_ci, dummy, opt_status)
      if (opt_status /= 0) status = 7

   end subroutine inv_transform_estimate

   subroutine minimize_transform_fit(x, y, family, jacobian_adjusted, weights, gamma, lower, upper, tolerance, &
                                     minimum, value, status)
      real(dp), intent(in) :: x(:, :) !! Design matrix passed to the transformation-fit objective.
      real(dp), intent(in) :: y(:) !! Response vector passed to the transformation-fit objective.
      character(len=*), intent(in) :: family !! Transformation family evaluated during the scalar search.
      logical, intent(in) :: jacobian_adjusted !! Whether objective evaluations apply Jacobian normalization.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights forwarded to the regression fit.
      real(dp), intent(in), optional :: gamma !! Optional BCN gamma parameter forwarded to the transformation.
      real(dp), intent(in) :: lower !! Lower bound of the lambda search interval.
      real(dp), intent(in) :: upper !! Upper bound of the lambda search interval; must exceed lower.
      real(dp), intent(in) :: tolerance !! Positive relative stopping tolerance for interval width.
      real(dp), intent(out) :: minimum !! Approximate lambda minimizing the transformation-fit objective.
      real(dp), intent(out) :: value !! Objective value at the returned lambda.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid bounds or nonfinite objective values.
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      integer :: iter

      if (upper <= lower .or. tolerance <= 0.0_dp) then
         minimum = ieee_value(0.0_dp, ieee_quiet_nan)
         value = minimum
         status = 1
         return
      end if
      a = lower
      b = upper
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = transform_fit_objective(x, y, family, c, jacobian_adjusted, weights, gamma)
      fd = transform_fit_objective(x, y, family, d, jacobian_adjusted, weights, gamma)
      if (.not. ieee_is_finite(fc) .or. .not. ieee_is_finite(fd)) then
         minimum = ieee_value(0.0_dp, ieee_quiet_nan)
         value = minimum
         status = 2
         return
      end if
      do iter = 1, 240
         if (abs(b - a) <= tolerance * max(1.0_dp, abs(a), abs(b))) exit
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = transform_fit_objective(x, y, family, c, jacobian_adjusted, weights, gamma)
            if (.not. ieee_is_finite(fc)) then
               minimum = ieee_value(0.0_dp, ieee_quiet_nan)
               value = minimum
               status = 2
               return
            end if
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = transform_fit_objective(x, y, family, d, jacobian_adjusted, weights, gamma)
            if (.not. ieee_is_finite(fd)) then
               minimum = ieee_value(0.0_dp, ieee_quiet_nan)
               value = minimum
               status = 2
               return
            end if
         end if
      end do
      minimum = 0.5_dp * (a + b)
      value = transform_fit_objective(x, y, family, minimum, jacobian_adjusted, weights, gamma)
      if (ieee_is_finite(value)) then
         status = 0
      else
         status = 2
      end if
   end subroutine minimize_transform_fit

   function inverse_rss_objective(x, y, family, lambda, gamma) result(value)
      real(dp), intent(in) :: x(:) !! Predictor values transformed for a straight-line inverse-response fit.
      real(dp), intent(in) :: y(:) !! Response values used in the straight-line fit.
      character(len=*), intent(in) :: family !! Transformation family applied to x.
      real(dp), intent(in) :: lambda !! Trial transformation power.
      real(dp), intent(in), optional :: gamma !! Optional BCN gamma parameter.
      real(dp) :: value
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: transformed(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: residuals(:)
      integer :: local_status

      allocate(design(size(x), 2), transformed(size(x)), beta(2), residuals(size(x)))
      call apply_transform(x, family, lambda, .false., transformed, local_status, gamma)
      if (local_status /= 0) then
         value = huge(1.0_dp)
         return
      end if
      design(:, 1) = 1.0_dp
      design(:, 2) = transformed
      call least_squares_residuals(design, y, beta, residuals, local_status)
      if (local_status /= 0) then
         value = huge(1.0_dp)
      else
         value = sum(residuals**2)
         if (value <= tiny(1.0_dp)) value = tiny(1.0_dp)
      end if
   end function inverse_rss_objective

   function inverse_ci_objective(x, y, family, lambda, dev0, sample_size, cutoff, gamma) result(value)
      real(dp), intent(in) :: x(:) !! Predictor values used by the inverse-response residual objective.
      real(dp), intent(in) :: y(:) !! Response values used by the inverse-response residual objective.
      character(len=*), intent(in) :: family !! Transformation family applied to x.
      real(dp), intent(in) :: lambda !! Trial transformation power for the confidence-limit equation.
      real(dp), intent(in) :: dev0 !! Deviance-like value at the fitted transformation power.
      real(dp), intent(in) :: sample_size !! Number of observations represented as real(dp).
      real(dp), intent(in) :: cutoff !! Half chi-square cutoff for the likelihood confidence boundary.
      real(dp), intent(in), optional :: gamma !! Optional BCN gamma parameter.
      real(dp) :: value
      real(dp) :: rss

      rss = inverse_rss_objective(x, y, family, lambda, gamma)
      if (.not. ieee_is_finite(rss) .or. rss >= huge(1.0_dp) / 2.0_dp) then
         value = huge(1.0_dp)
      else
         value = abs(dev0 + sample_size * log(rss) - cutoff)
      end if
   end function inverse_ci_objective

   subroutine minimize_inverse_rss(x, y, family, gamma, lower, upper, tolerance, minimum, value, status)
      real(dp), intent(in) :: x(:) !! Predictor values used by the inverse-response residual objective.
      real(dp), intent(in) :: y(:) !! Response values used by the inverse-response residual objective.
      character(len=*), intent(in) :: family !! Transformation family applied during the scalar search.
      real(dp), intent(in), optional :: gamma !! Optional BCN gamma parameter.
      real(dp), intent(in) :: lower !! Lower bound of the lambda search interval.
      real(dp), intent(in) :: upper !! Upper bound of the lambda search interval; must exceed lower.
      real(dp), intent(in) :: tolerance !! Positive relative stopping tolerance for interval width.
      real(dp), intent(out) :: minimum !! Approximate lambda minimizing residual sum of squares.
      real(dp), intent(out) :: value !! Residual sum of squares at the returned lambda.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid bounds or nonfinite objective values.
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      integer :: iter

      if (upper <= lower .or. tolerance <= 0.0_dp) then
         minimum = ieee_value(0.0_dp, ieee_quiet_nan)
         value = minimum
         status = 1
         return
      end if
      a = lower
      b = upper
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = inverse_rss_objective(x, y, family, c, gamma)
      fd = inverse_rss_objective(x, y, family, d, gamma)
      if (.not. ieee_is_finite(fc) .or. .not. ieee_is_finite(fd)) then
         minimum = ieee_value(0.0_dp, ieee_quiet_nan)
         value = minimum
         status = 2
         return
      end if
      do iter = 1, 240
         if (abs(b - a) <= tolerance * max(1.0_dp, abs(a), abs(b))) exit
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = inverse_rss_objective(x, y, family, c, gamma)
            if (.not. ieee_is_finite(fc)) then
               minimum = ieee_value(0.0_dp, ieee_quiet_nan)
               value = minimum
               status = 2
               return
            end if
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = inverse_rss_objective(x, y, family, d, gamma)
            if (.not. ieee_is_finite(fd)) then
               minimum = ieee_value(0.0_dp, ieee_quiet_nan)
               value = minimum
               status = 2
               return
            end if
         end if
      end do
      minimum = 0.5_dp * (a + b)
      value = inverse_rss_objective(x, y, family, minimum, gamma)
      if (ieee_is_finite(value)) then
         status = 0
      else
         status = 2
      end if
   end subroutine minimize_inverse_rss

   subroutine minimize_inverse_ci(x, y, family, gamma, dev0, sample_size, cutoff, lower, upper, tolerance, &
                                  minimum, value, status)
      real(dp), intent(in) :: x(:) !! Predictor values used by the inverse-response confidence objective.
      real(dp), intent(in) :: y(:) !! Response values used by the inverse-response confidence objective.
      character(len=*), intent(in) :: family !! Transformation family applied during the scalar search.
      real(dp), intent(in), optional :: gamma !! Optional BCN gamma parameter.
      real(dp), intent(in) :: dev0 !! Deviance-like value at the fitted transformation power.
      real(dp), intent(in) :: sample_size !! Number of observations represented as real(dp).
      real(dp), intent(in) :: cutoff !! Half chi-square cutoff defining the likelihood confidence boundary.
      real(dp), intent(in) :: lower !! Lower bound of the lambda search interval.
      real(dp), intent(in) :: upper !! Upper bound of the lambda search interval; must exceed lower.
      real(dp), intent(in) :: tolerance !! Positive relative stopping tolerance for interval width.
      real(dp), intent(out) :: minimum !! Approximate lambda minimizing the confidence-boundary objective.
      real(dp), intent(out) :: value !! Absolute confidence-boundary equation residual at the returned lambda.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid bounds or nonfinite objective values.
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      integer :: iter

      if (upper <= lower .or. tolerance <= 0.0_dp) then
         minimum = ieee_value(0.0_dp, ieee_quiet_nan)
         value = minimum
         status = 1
         return
      end if
      a = lower
      b = upper
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = inverse_ci_objective(x, y, family, c, dev0, sample_size, cutoff, gamma)
      fd = inverse_ci_objective(x, y, family, d, dev0, sample_size, cutoff, gamma)
      if (.not. ieee_is_finite(fc) .or. .not. ieee_is_finite(fd)) then
         minimum = ieee_value(0.0_dp, ieee_quiet_nan)
         value = minimum
         status = 2
         return
      end if
      do iter = 1, 240
         if (abs(b - a) <= tolerance * max(1.0_dp, abs(a), abs(b))) exit
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = inverse_ci_objective(x, y, family, c, dev0, sample_size, cutoff, gamma)
            if (.not. ieee_is_finite(fc)) then
               minimum = ieee_value(0.0_dp, ieee_quiet_nan)
               value = minimum
               status = 2
               return
            end if
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = inverse_ci_objective(x, y, family, d, dev0, sample_size, cutoff, gamma)
            if (.not. ieee_is_finite(fd)) then
               minimum = ieee_value(0.0_dp, ieee_quiet_nan)
               value = minimum
               status = 2
               return
            end if
         end if
      end do
      minimum = 0.5_dp * (a + b)
      value = inverse_ci_objective(x, y, family, minimum, dev0, sample_size, cutoff, gamma)
      if (ieee_is_finite(value)) then
         status = 0
      else
         status = 2
      end if
   end subroutine minimize_inverse_ci

   function transform_fit_objective(x, y, family, lambda, jacobian_adjusted, weights, gamma) result(value)
      real(dp), intent(in) :: x(:, :) !! Design matrix used to residualize the transformed response.
      real(dp), intent(in) :: y(:) !! Response vector to transform.
      character(len=*), intent(in) :: family !! Transformation family name.
      real(dp), intent(in) :: lambda !! Trial transformation parameter.
      logical, intent(in) :: jacobian_adjusted !! Whether to apply the family Jacobian normalization.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights.
      real(dp), intent(in), optional :: gamma !! Optional BCN gamma parameter.
      real(dp) :: value
      real(dp), allocatable :: transformed(:), beta(:), residuals(:)
      real(dp) :: variance
      integer :: transform_status, fit_status, variance_status, n_valid, n

      n = size(y)
      allocate(transformed(n), beta(size(x, 2)), residuals(n))
      call apply_transform(y, family, lambda, jacobian_adjusted, transformed, transform_status, gamma)
      if (transform_status /= 0) then
         value = huge(1.0_dp)
         return
      end if
      if (present(weights)) then
         call least_squares_residuals(x, transformed, beta, residuals, fit_status, weights)
      else
         call least_squares_residuals(x, transformed, beta, residuals, fit_status)
      end if
      if (fit_status /= 0) then
         value = huge(1.0_dp)
         return
      end if
      call sample_variance(residuals, variance, n_valid, variance_status)
      if (variance_status /= 0 .or. variance <= tiny(1.0_dp)) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.5_dp * real(n, dp) * log(real(n - 1, dp) / real(n, dp) * variance)
   end function transform_fit_objective

   subroutine apply_transform(values, family, lambda, jacobian_adjusted, transformed, status, gamma)
      real(dp), intent(in) :: values(:) !! Values transformed using the selected power family.
      character(len=*), intent(in) :: family !! Family name matched case-insensitively: bcPower, yjPower, or bcnPower.
      real(dp), intent(in) :: lambda !! Power parameter for the selected family.
      logical, intent(in) :: jacobian_adjusted !! Whether to apply the family Jacobian normalization.
      real(dp), intent(out) :: transformed(:) !! Transformed values with the same length as values.
      integer, intent(out) :: status !! Zero on success; nonzero for unknown family or transformation failure.
      real(dp), intent(in), optional :: gamma !! BCN gamma parameter; omitted gamma means zero.
      character(len=16) :: kind
      real(dp) :: gamma_value

      kind = lowercase(trim(family))
      select case (trim(kind))
      case ("bcpower", "boxcox", "bc")
         call bc_power(values, lambda, transformed, status, jacobian_adjusted)
      case ("yjpower", "yeojohnson", "yj")
         call yj_power(values, lambda, transformed, status, jacobian_adjusted)
      case ("bcnpower", "bcn")
         gamma_value = 0.0_dp
         if (present(gamma)) gamma_value = gamma
         call bcn_power(values, lambda, gamma_value, transformed, status, jacobian_adjusted)
      case default
         transformed = ieee_value(0.0_dp, ieee_quiet_nan)
         status = 9
      end select
   end subroutine apply_transform

   integer function validate_fit_input(x, y, weights) result(status)
      real(dp), intent(in) :: x(:, :) !! Design matrix to validate for a transformation fit.
      real(dp), intent(in) :: y(:) !! Response vector to validate for a transformation fit.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights.

      status = 0
      if (size(x, 1) /= size(y) .or. size(y) < 2 .or. size(x, 2) < 1) then
         status = 1
         return
      end if
      if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) then
         status = 2
         return
      end if
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp) .or. any(ieee_is_nan(weights))) status = 3
      end if
   end function validate_fit_input

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text !! Text converted to lowercase for family matching.
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

end module car_power_fit
