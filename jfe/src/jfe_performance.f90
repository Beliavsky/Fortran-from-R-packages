! SPDX-License-Identifier: GPL-2.0-or-later
module jfe_performance
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jfe_kinds
   use jfe_stats
   implicit none
   private

   type, public :: annualized_summary
      real(dp) :: annualized_return = 0.0_dp
      real(dp) :: annualized_sd = 0.0_dp
      real(dp) :: annualized_sharpe = 0.0_dp
      integer :: status = jfe_success
   end type annualized_summary

   type, public :: durbin_h_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: durbin_watson = 0.0_dp
      integer :: status = jfe_success
   end type durbin_h_result

   public :: active_premium, adjusted_sharpe_ratio, bernardo_ledoit_ratio
   public :: burke_ratio, d_ratio, kelly_ratio, martin_ratio
   public :: skewness_kurtosis_ratio, pain_index, mean_absolute_deviation
   public :: calmar_ratio, sterling_ratio, appraisal_ratio, tracking_error
   public :: information_ratio, treynor_ratio, downside_deviation
   public :: downside_potential, upside_risk, omega_sharpe_ratio, sortino_ratio
   public :: prospect_ratio, volatility_skewness, m2_sortino
   public :: sharpe_ratio, sharpe_ratio_annualized, pain_ratio
   public :: table_annualized_returns, return_annualized, capm_jensen_alpha
   public :: ulcer_index, drawdown_peak, max_drawdown, drawdowns
   public :: skewness, kurtosis, value_at_risk, expected_shortfall
   public :: durbin_h

contains

   subroutine set_status(optional_status, value)
      integer, intent(out), optional :: optional_status
      integer, intent(in) :: value
      if (present(optional_status)) optional_status = value
   end subroutine set_status

   real(dp) function return_annualized(r, scale, geometric, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      logical, intent(in), optional :: geometric
      integer, intent(out), optional :: status
      integer :: istat

      ans = annualized_return_value(r, scale, geometric, istat)
      call set_status(status, istat)
   end function return_annualized

   real(dp) function active_premium(ra, rb, scale, geometric, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      logical, intent(in), optional :: geometric
      integer, intent(out), optional :: status
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: aa, ab
      integer :: s1, s2

      call clean_pair(ra, rb, a, b)
      aa = annualized_return_value(a, scale, geometric, s1)
      ab = annualized_return_value(b, scale, geometric, s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         ans = aa - ab
         call set_status(status, jfe_success)
      end if
   end function active_premium

   real(dp) function value_at_risk(r, alpha, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: alpha
      integer, intent(out), optional :: status
      real(dp) :: a
      integer :: istat

      a = 0.05_dp
      if (present(alpha)) a = alpha
      ans = quantile_type1(r, a, istat)
      call set_status(status, istat)
   end function value_at_risk

   real(dp) function expected_shortfall(r, alpha, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: alpha
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: a, q
      integer :: istat

      a = 0.05_dp
      if (present(alpha)) a = alpha
      if (a <= 0.0_dp .or. a > 1.0_dp) then
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end if
      call clean_vector(r, x)
      q = quantile_type1(x, a, istat)
      if (istat /= jfe_success) then
         ans = nan_dp()
         call set_status(status, istat)
         return
      end if
      ! Algebraically identical to the package's CVaR expression.
      ans = q - sum(max(q - x, 0.0_dp))/(a*real(size(x), dp))
      call set_status(status, jfe_success)
   end function expected_shortfall

   real(dp) function sharpe_ratio(r, rf, alpha, risk_method, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: rf, alpha
      integer, intent(in), optional :: risk_method
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: risk_free, a, numerator, denominator
      integer :: method, istat

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      a = 0.05_dp
      if (present(alpha)) a = alpha
      method = risk_stddev
      if (present(risk_method)) method = risk_method
      call clean_vector(r, x)
      if (size(x) < 2) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
         return
      end if
      numerator = sum(x - risk_free)/real(size(x), dp)
      select case (method)
      case (risk_stddev)
         denominator = sd_value(x, istat)
      case (risk_var)
         denominator = -value_at_risk(x - risk_free, a, istat)
      case (risk_es)
         denominator = -expected_shortfall(x - risk_free, a, istat)
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      if (istat /= jfe_success .or. abs(denominator) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = numerator/denominator
         call set_status(status, jfe_success)
      end if
   end function sharpe_ratio

   real(dp) function sharpe_ratio_annualized(r, rf, alpha, scale, geometric, &
      risk_method, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: rf, alpha
      logical, intent(in), optional :: geometric
      integer, intent(in), optional :: risk_method
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: risk_free, a, numerator, denominator
      integer :: method, istat

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      a = 0.05_dp
      if (present(alpha)) a = alpha
      method = risk_stddev
      if (present(risk_method)) method = risk_method
      call clean_vector(r, x)
      numerator = annualized_return_value(x, scale, geometric, istat) - risk_free
      if (istat /= jfe_success) then
         ans = nan_dp()
         call set_status(status, istat)
         return
      end if
      select case (method)
      case (risk_stddev)
         denominator = sd_value(x, istat)*sqrt(scale)
      case (risk_var)
         ! The upstream implementation does not annualize historical VaR.
         denominator = -value_at_risk(x, a, istat)
      case (risk_es)
         ! The upstream implementation does not annualize historical ES.
         denominator = -expected_shortfall(x, a, istat)
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      if (istat /= jfe_success .or. abs(denominator) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = numerator/denominator
         call set_status(status, jfe_success)
      end if
   end function sharpe_ratio_annualized

   real(dp) function adjusted_sharpe_ratio(r, rf, scale, risk_method, alpha, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: rf, alpha
      integer, intent(in), optional :: risk_method
      integer, intent(out), optional :: status
      real(dp) :: sr, sk, ku
      integer :: s1, s2, s3

      sr = sharpe_ratio_annualized(r, rf, alpha, scale, risk_method=risk_method, status=s1)
      sk = skewness_value(r, skew_moment, s2)
      ku = kurtosis_value(r, kurt_excess, s3)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. s3 /= jfe_success) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         ans = sr*(1.0_dp + sr*sk/6.0_dp - sr*sr*ku/24.0_dp)
         call set_status(status, jfe_success)
      end if
   end function adjusted_sharpe_ratio

   real(dp) function bernardo_ledoit_ratio(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: losses

      call clean_vector(r, x)
      losses = -sum(x, mask=x < 0.0_dp)
      if (size(x) < 1) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else if (losses <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = sum(x, mask=x > 0.0_dp)/losses
         call set_status(status, jfe_success)
      end if
   end function bernardo_ledoit_ratio

   real(dp) function burke_ratio(r, rf, scale, modified, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: rf
      logical, intent(in), optional :: modified
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:), episodes(:)
      real(dp) :: risk_free, temp, denominator, annual_return
      integer :: i, j, n_episode, start, istat
      logical :: in_drawdown, use_modified

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      use_modified = .false.
      if (present(modified)) use_modified = modified
      call clean_vector(r, x)
      if (size(x) < 1) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
         return
      end if
      allocate(episodes(size(x)))
      episodes = 0.0_dp
      n_episode = 0
      in_drawdown = .false.
      start = 1
      do i = 2, size(x)
         if (x(i) < 0.0_dp) then
            if (.not. in_drawdown) then
               start = i
               in_drawdown = .true.
            end if
         else if (in_drawdown) then
            temp = 1.0_dp
            do j = start, i - 1
               temp = temp*(1.0_dp + x(j)*0.01_dp)
            end do
            n_episode = n_episode + 1
            episodes(n_episode) = (temp - 1.0_dp)*100.0_dp
            in_drawdown = .false.
         end if
      end do
      if (in_drawdown) then
         temp = 1.0_dp
         do j = start, size(x)
            temp = temp*(1.0_dp + x(j)*0.01_dp)
         end do
         n_episode = n_episode + 1
         episodes(n_episode) = (temp - 1.0_dp)*100.0_dp
      end if
      if (n_episode == 0) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
         return
      end if
      denominator = sqrt(sum(episodes(1:n_episode)**2))
      if (use_modified) denominator = denominator/sqrt(real(size(x), dp))
      annual_return = annualized_return_value(x, scale, .true., istat)
      if (istat /= jfe_success .or. denominator <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = (annual_return - risk_free)/denominator
         call set_status(status, jfe_success)
      end if
   end function burke_ratio

   real(dp) function d_ratio(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: den
      integer :: nd, nu

      call clean_vector(r, x)
      nd = count(x < 0.0_dp)
      nu = count(x > 0.0_dp)
      den = real(nu, dp)*sum(x, mask=x > 0.0_dp)
      if (nd == 0 .or. nu == 0 .or. abs(den) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = -real(nd, dp)*sum(x, mask=x < 0.0_dp)/den
         call set_status(status, jfe_success)
      end if
   end function d_ratio

   real(dp) function kelly_ratio(r, rf, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: rf
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: risk_free, varr
      integer :: istat

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      call clean_vector(r, x)
      varr = variance_value(x, istat)
      if (istat /= jfe_success .or. varr <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = (sum(x - risk_free)/real(size(x), dp))/varr/2.0_dp
         call set_status(status, jfe_success)
      end if
   end function kelly_ratio

   function drawdown_peak(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      integer :: istat

      ans = drawdown_peak_value(r, istat)
      call set_status(status, istat)
   end function drawdown_peak

   function drawdowns(r, geometric, status) result(ans)
      real(dp), intent(in) :: r(:)
      logical, intent(in), optional :: geometric
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      integer :: istat

      ans = drawdowns_value(r, geometric, istat)
      call set_status(status, istat)
   end function drawdowns

   real(dp) function max_drawdown(r, geometric, invert, status) result(ans)
      real(dp), intent(in) :: r(:)
      logical, intent(in), optional :: geometric, invert
      integer, intent(out), optional :: status
      integer :: istat

      ans = max_drawdown_value(r, geometric, invert, istat)
      call set_status(status, istat)
   end function max_drawdown

   real(dp) function ulcer_index(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: dd(:), x(:)

      dd = drawdown_peak_value(r)
      call clean_vector(dd, x)
      if (size(x) < 1) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         ans = sqrt(sum(x*x)/real(size(x), dp))
         call set_status(status, jfe_success)
      end if
   end function ulcer_index

   real(dp) function pain_index(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: dd(:), x(:)

      dd = drawdown_peak_value(r)
      call clean_vector(dd, x)
      if (size(x) < 1) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         ans = sum(abs(x))/real(size(x), dp)
         call set_status(status, jfe_success)
      end if
   end function pain_index

   real(dp) function martin_ratio(r, rf, scale, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: rf
      integer, intent(out), optional :: status
      real(dp) :: risk_free, ui, ar
      integer :: s1, s2

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      ui = ulcer_index(r, s1)
      ar = annualized_return_value(r, scale, .true., s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(ui) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = (ar - risk_free)/ui
         call set_status(status, jfe_success)
      end if
   end function martin_ratio

   real(dp) function skewness_kurtosis_ratio(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp) :: sk, ku
      integer :: s1, s2

      sk = skewness_value(r, skew_moment, s1)
      ku = kurtosis_value(r, kurt_moment, s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(ku) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = sk/ku
         call set_status(status, jfe_success)
      end if
   end function skewness_kurtosis_ratio

   real(dp) function mean_absolute_deviation(r, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: mu

      call clean_vector(r, x)
      if (size(x) < 1) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         mu = sum(x)/real(size(x), dp)
         ans = sum(abs(x - mu))/real(size(x), dp)
         call set_status(status, jfe_success)
      end if
   end function mean_absolute_deviation

   real(dp) function calmar_ratio(r, scale, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      integer, intent(out), optional :: status
      real(dp) :: ar, mdd
      integer :: s1, s2

      ar = annualized_return_value(r, scale, .true., s1)
      mdd = max_drawdown_value(r, .true., .true., s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(mdd) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = ar/abs(mdd)
         call set_status(status, jfe_success)
      end if
   end function calmar_ratio

   real(dp) function sterling_ratio(r, scale, excess, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: excess
      integer, intent(out), optional :: status
      real(dp) :: ex, ar, den
      integer :: s1, s2

      ex = 0.1_dp
      if (present(excess)) ex = excess
      ar = annualized_return_value(r, scale, .true., s1)
      den = abs(max_drawdown_value(r, .true., .false., s2) + ex)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. den <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = ar/den
         call set_status(status, jfe_success)
      end if
   end function sterling_ratio

   real(dp) function capm_jensen_alpha(ra, rb, rf, scale, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      real(dp), intent(in), optional :: rf
      integer, intent(out), optional :: status
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: risk_free, beta, ar_a, ar_b
      integer :: s1, s2, s3

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      call clean_pair(ra, rb, a, b)
      beta = beta_value(a - risk_free, b - risk_free, s1)
      ar_a = annualized_return_value(a, scale, .true., s2)
      ar_b = annualized_return_value(b, scale, .true., s3)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. s3 /= jfe_success) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         ans = ar_a - risk_free - beta*(ar_b - risk_free)
         call set_status(status, jfe_success)
      end if
   end function capm_jensen_alpha

   real(dp) function appraisal_ratio(ra, rb, rf, scale, method, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      real(dp), intent(in), optional :: rf
      integer, intent(in), optional :: method
      integer, intent(out), optional :: status
      real(dp), allocatable :: a(:), b(:), residuals(:)
      real(dp) :: risk_free, alpha, beta, intercept, specific_risk, systematic_risk
      integer :: meth, s1, s2

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      meth = appraisal_standard
      if (present(method)) meth = method
      call clean_pair(ra, rb, a, b)
      alpha = capm_jensen_alpha(a, b, risk_free, scale, s1)
      call ols_fit(a - risk_free, b - risk_free, intercept, beta, residuals, s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
         return
      end if
      select case (meth)
      case (appraisal_standard)
         specific_risk = sqrt(sum((residuals - sum(residuals)/real(size(residuals), dp))**2)/ &
            real(size(residuals), dp))*sqrt(scale)
         if (specific_risk <= tiny(1.0_dp)) then
            ans = nan_dp()
            call set_status(status, jfe_zero_denominator)
            return
         end if
         ans = alpha/specific_risk
      case (appraisal_modified)
         if (abs(beta) <= tiny(1.0_dp)) then
            ans = nan_dp()
            call set_status(status, jfe_zero_denominator)
            return
         end if
         ans = alpha/beta
      case (appraisal_alternative)
         systematic_risk = beta*sd_value(b)*sqrt(scale)
         if (abs(systematic_risk) <= tiny(1.0_dp)) then
            ans = nan_dp()
            call set_status(status, jfe_zero_denominator)
            return
         end if
         ans = alpha/systematic_risk
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      call set_status(status, jfe_success)
   end function appraisal_ratio

   real(dp) function tracking_error(ra, rb, scale, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      integer, intent(out), optional :: status
      real(dp), allocatable :: a(:), b(:)
      integer :: istat

      call clean_pair(ra, rb, a, b)
      ans = sd_value(a - b, istat)*sqrt(scale)
      call set_status(status, istat)
   end function tracking_error

   real(dp) function information_ratio(ra, rb, scale, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      integer, intent(out), optional :: status
      real(dp) :: ap, te
      integer :: s1, s2

      ap = active_premium(ra, rb, scale, status=s1)
      te = tracking_error(ra, rb, scale, s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(te) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = ap/te
         call set_status(status, jfe_success)
      end if
   end function information_ratio

   real(dp) function treynor_ratio(ra, rb, rf, scale, modified, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      real(dp), intent(in), optional :: rf
      logical, intent(in), optional :: modified
      integer, intent(out), optional :: status
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: risk_free, beta, numerator, den
      integer :: s1, s2
      logical :: use_modified

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      use_modified = .false.
      if (present(modified)) use_modified = modified
      call clean_pair(ra, rb, a, b)
      beta = beta_value(a - risk_free, b - risk_free, s1)
      numerator = annualized_return_value(a - risk_free, scale, .true., s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
         return
      end if
      if (use_modified) then
         den = beta*sd_value(b - risk_free)*sqrt(scale)
      else
         den = beta
      end if
      if (abs(den) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = numerator/den
         call set_status(status, jfe_success)
      end if
   end function treynor_ratio

   real(dp) function downside_deviation(r, mar, method, potential, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: mar
      integer, intent(in), optional :: method
      logical, intent(in), optional :: potential
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:), below(:)
      real(dp) :: threshold, divisor
      integer :: meth, i, n
      logical :: use_potential

      threshold = 0.0_dp
      if (present(mar)) threshold = mar
      meth = downside_full
      if (present(method)) meth = method
      use_potential = .false.
      if (present(potential)) use_potential = potential
      call clean_vector(r, x)
      n = count(x < threshold)
      allocate(below(n))
      n = 0
      do i = 1, size(x)
         if (x(i) < threshold) then
            n = n + 1
            below(n) = x(i)
         end if
      end do
      if (size(x) < 1) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
         return
      end if
      select case (meth)
      case (downside_full)
         divisor = real(size(x), dp)
      case (downside_subset)
         if (size(below) < 1) then
            ans = 0.0_dp
            call set_status(status, jfe_success)
            return
         end if
         divisor = real(size(below), dp)
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      if (use_potential) then
         ans = sum(threshold - below)/divisor
      else
         ans = sqrt(sum((threshold - below)**2)/divisor)
      end if
      call set_status(status, jfe_success)
   end function downside_deviation

   real(dp) function downside_potential(r, mar, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: mar
      integer, intent(out), optional :: status
      integer :: istat

      ans = downside_deviation(r, mar, downside_full, .true., istat)
      call set_status(status, istat)
   end function downside_potential

   real(dp) function upside_risk(r, mar, method, statistic, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: mar
      integer, intent(in), optional :: method, statistic
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:), above(:)
      real(dp) :: threshold, divisor, second_moment
      integer :: meth, stat, i, n

      threshold = 0.0_dp
      if (present(mar)) threshold = mar
      meth = downside_full
      if (present(method)) meth = method
      stat = 1
      if (present(statistic)) stat = statistic
      call clean_vector(r, x)
      n = count(x > threshold)
      allocate(above(n))
      n = 0
      do i = 1, size(x)
         if (x(i) > threshold) then
            n = n + 1
            above(n) = x(i)
         end if
      end do
      select case (meth)
      case (downside_full)
         divisor = real(size(x), dp)
      case (downside_subset)
         if (size(above) == 0) then
            ans = 0.0_dp
            call set_status(status, jfe_success)
            return
         end if
         divisor = real(size(above), dp)
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      if (divisor <= 0.0_dp) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
         return
      end if
      second_moment = sum((above - threshold)**2)/divisor
      select case (stat)
      case (1)
         ans = sqrt(second_moment)
      case (2)
         ans = second_moment
      case (3)
         ans = sum(above - threshold)/divisor
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      call set_status(status, jfe_success)
   end function upside_risk

   real(dp) function omega_sharpe_ratio(r, mar, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: mar
      integer, intent(out), optional :: status
      real(dp) :: up, down
      integer :: s1, s2

      up = upside_risk(r, mar, downside_full, 3, s1)
      down = downside_potential(r, mar, s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(down) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = (up - down)/down
         call set_status(status, jfe_success)
      end if
   end function omega_sharpe_ratio

   real(dp) function sortino_ratio(r, mar, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: mar
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: threshold, dd
      integer :: istat

      threshold = 0.0_dp
      if (present(mar)) threshold = mar
      call clean_vector(r, x)
      dd = downside_deviation(x, threshold, status=istat)
      if (istat /= jfe_success .or. abs(dd) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = sum(x - threshold)/real(size(x), dp)/dd
         call set_status(status, jfe_success)
      end if
   end function sortino_ratio

   real(dp) function prospect_ratio(r, mar, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in) :: mar
      integer, intent(out), optional :: status
      real(dp), allocatable :: x(:)
      real(dp) :: dd, numerator
      integer :: istat

      call clean_vector(r, x)
      dd = downside_deviation(x, mar, status=istat)
      if (istat /= jfe_success .or. abs(dd) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         numerator = sum(x, mask=x > 0.0_dp) + 2.25_dp*sum(x, mask=x < 0.0_dp) - mar
         ans = numerator/(dd*real(size(x), dp))
         call set_status(status, jfe_success)
      end if
   end function prospect_ratio

   real(dp) function volatility_skewness(r, mar, statistic, status) result(ans)
      real(dp), intent(in) :: r(:)
      real(dp), intent(in), optional :: mar
      integer, intent(in), optional :: statistic
      integer, intent(out), optional :: status
      real(dp) :: up, down
      integer :: stat, s1, s2

      stat = volatility_ratio
      if (present(statistic)) stat = statistic
      select case (stat)
      case (volatility_ratio)
         up = upside_risk(r, mar, downside_full, 2, s1)
         down = downside_deviation(r, mar, status=s2)**2
      case (variability_ratio)
         up = upside_risk(r, mar, downside_full, 1, s1)
         down = downside_deviation(r, mar, status=s2)
      case default
         ans = nan_dp()
         call set_status(status, jfe_invalid_argument)
         return
      end select
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(down) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = up/down
         call set_status(status, jfe_success)
      end if
   end function volatility_skewness

   real(dp) function m2_sortino(ra, rb, mar, scale, status) result(ans)
      real(dp), intent(in) :: ra(:), rb(:), scale
      real(dp), intent(in), optional :: mar
      integer, intent(out), optional :: status
      real(dp) :: threshold, rp, sigmad, sigmadm, sr
      integer :: s1, s2, s3, s4

      threshold = 0.0_dp
      if (present(mar)) threshold = mar
      rp = annualized_return_value(ra, scale, .true., s1)
      sigmad = downside_deviation(ra, threshold, status=s2)*sqrt(scale)
      sigmadm = downside_deviation(rb, threshold, status=s3)*sqrt(scale)
      sr = sortino_ratio(ra, threshold, s4)
      if (any([s1, s2, s3, s4] /= jfe_success)) then
         ans = nan_dp()
         call set_status(status, jfe_insufficient_data)
      else
         ans = rp + sr*(sigmadm - sigmad)
         call set_status(status, jfe_success)
      end if
   end function m2_sortino

   real(dp) function pain_ratio(r, rf, scale, status) result(ans)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: rf
      integer, intent(out), optional :: status
      real(dp) :: risk_free, pi, ar
      integer :: s1, s2

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      pi = pain_index(r, s1)
      ar = annualized_return_value(r, scale, .true., s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success .or. abs(pi) <= tiny(1.0_dp)) then
         ans = nan_dp()
         call set_status(status, jfe_zero_denominator)
      else
         ans = (ar - risk_free)/pi
         call set_status(status, jfe_success)
      end if
   end function pain_ratio

   function table_annualized_returns(r, scale, rf, geometric) result(summary)
      real(dp), intent(in) :: r(:), scale
      real(dp), intent(in), optional :: rf
      logical, intent(in), optional :: geometric
      type(annualized_summary) :: summary
      real(dp) :: risk_free
      integer :: s1, s2, s3

      risk_free = 0.0_dp
      if (present(rf)) risk_free = rf
      summary%annualized_return = annualized_return_value(r, scale, geometric, s1)
      summary%annualized_sd = sd_value(r, s2)*sqrt(scale)
      summary%annualized_sharpe = sharpe_ratio_annualized(r, risk_free, scale=scale, &
         geometric=geometric, status=s3)
      if (any([s1, s2, s3] /= jfe_success)) then
         summary%status = jfe_insufficient_data
      else
         summary%status = jfe_success
      end if
   end function table_annualized_returns

   real(dp) function skewness(r, method, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(in), optional :: method
      integer, intent(out), optional :: status
      integer :: istat

      ans = skewness_value(r, method, istat)
      call set_status(status, istat)
   end function skewness

   real(dp) function kurtosis(r, method, status) result(ans)
      real(dp), intent(in) :: r(:)
      integer, intent(in), optional :: method
      integer, intent(out), optional :: status
      integer :: istat

      ans = kurtosis_value(r, method, istat)
      call set_status(status, istat)
   end function kurtosis

   function durbin_h(residuals, n_fitted, n_coefficients, lag_variance) result(ans)
      real(dp), intent(in) :: residuals(:)
      integer, intent(in) :: n_fitted, n_coefficients
      real(dp), intent(in) :: lag_variance
      type(durbin_h_result) :: ans
      real(dp), allocatable :: e(:)
      real(dp) :: denominator, hden
      integer :: n

      call clean_vector(residuals, e)
      if (size(e) < 2 .or. n_fitted < 1 .or. n_coefficients < 1 .or. lag_variance < 0.0_dp) then
         ans%statistic = nan_dp()
         ans%p_value = nan_dp()
         ans%durbin_watson = nan_dp()
         ans%status = jfe_invalid_argument
         return
      end if
      denominator = sum(e*e)
      if (denominator <= tiny(1.0_dp)) then
         ans%statistic = nan_dp()
         ans%p_value = nan_dp()
         ans%durbin_watson = nan_dp()
         ans%status = jfe_zero_denominator
         return
      end if
      ans%durbin_watson = sum((e(2:) - e(:size(e) - 1))**2)/denominator
      n = n_fitted - n_coefficients + 1
      hden = 1.0_dp - real(n, dp)*lag_variance
      if (n <= 0 .or. hden <= 0.0_dp) then
         ans%statistic = nan_dp()
         ans%p_value = nan_dp()
         ans%status = jfe_invalid_argument
         return
      end if
      ans%statistic = (1.0_dp - ans%durbin_watson/2.0_dp)*sqrt(real(n, dp)/hden)
      ans%p_value = 2.0_dp*(1.0_dp - normal_cdf(abs(ans%statistic)))
      ans%status = jfe_success
   end function durbin_h

end module jfe_performance
