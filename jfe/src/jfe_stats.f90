! SPDX-License-Identifier: GPL-2.0-or-later
module jfe_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use jfe_kinds, only : dp, jfe_success, jfe_invalid_argument, jfe_insufficient_data, &
      jfe_zero_denominator, skew_moment, skew_fisher, skew_sample, kurt_excess, &
      kurt_moment, kurt_fisher, kurt_sample, kurt_sample_excess
   implicit none
   private

   public :: clean_vector, clean_pair, nan_dp
   public :: mean_value, variance_value, sd_value, covariance_value, beta_value
   public :: quantile_type1, skewness_value, kurtosis_value, normal_cdf
   public :: annualized_return_value, drawdowns_value, drawdown_peak_value
   public :: max_drawdown_value, ols_fit

contains

   pure real(dp) function nan_dp() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   subroutine clean_vector(x, y)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: y(:)
      integer :: i, n

      n = count(ieee_is_finite(x))
      allocate(y(n))
      n = 0
      do i = 1, size(x)
         if (ieee_is_finite(x(i))) then
            n = n + 1
            y(n) = x(i)
         end if
      end do
   end subroutine clean_vector

   subroutine clean_pair(x, y, xc, yc)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), allocatable, intent(out) :: xc(:), yc(:)
      integer :: i, n, m

      m = min(size(x), size(y))
      n = 0
      do i = 1, m
         if (ieee_is_finite(x(i)) .and. ieee_is_finite(y(i))) n = n + 1
      end do
      allocate(xc(n), yc(n))
      n = 0
      do i = 1, m
         if (ieee_is_finite(x(i)) .and. ieee_is_finite(y(i))) then
            n = n + 1
            xc(n) = x(i)
            yc(n) = y(i)
         end if
      end do
   end subroutine clean_pair

   real(dp) function mean_value(x, status) result(mu)
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)

      call clean_vector(x, y)
      if (size(y) < 1) then
         mu = nan_dp()
         if (present(status)) status = jfe_insufficient_data
      else
         mu = sum(y)/real(size(y), dp)
         if (present(status)) status = jfe_success
      end if
   end function mean_value

   real(dp) function variance_value(x, status) result(v)
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      real(dp) :: mu

      call clean_vector(x, y)
      if (size(y) < 2) then
         v = nan_dp()
         if (present(status)) status = jfe_insufficient_data
         return
      end if
      mu = sum(y)/real(size(y), dp)
      v = sum((y - mu)**2)/real(size(y) - 1, dp)
      if (present(status)) status = jfe_success
   end function variance_value

   real(dp) function sd_value(x, status) result(s)
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      integer :: istat
      real(dp) :: v

      v = variance_value(x, istat)
      if (istat /= jfe_success .or. v < 0.0_dp) then
         s = nan_dp()
      else
         s = sqrt(v)
      end if
      if (present(status)) status = istat
   end function sd_value

   real(dp) function covariance_value(x, y, status) result(covxy)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: xc(:), yc(:)
      real(dp) :: mx, my

      call clean_pair(x, y, xc, yc)
      if (size(xc) < 2) then
         covxy = nan_dp()
         if (present(status)) status = jfe_insufficient_data
         return
      end if
      mx = sum(xc)/real(size(xc), dp)
      my = sum(yc)/real(size(yc), dp)
      covxy = sum((xc - mx)*(yc - my))/real(size(xc) - 1, dp)
      if (present(status)) status = jfe_success
   end function covariance_value

   real(dp) function beta_value(y, x, status) result(beta)
      real(dp), intent(in) :: y(:), x(:)
      integer, intent(out), optional :: status
      real(dp) :: vx, cv
      integer :: s1, s2

      vx = variance_value(x, s1)
      cv = covariance_value(y, x, s2)
      if (s1 /= jfe_success .or. s2 /= jfe_success) then
         beta = nan_dp()
         if (present(status)) status = jfe_insufficient_data
      else if (abs(vx) <= tiny(1.0_dp)) then
         beta = nan_dp()
         if (present(status)) status = jfe_zero_denominator
      else
         beta = cv/vx
         if (present(status)) status = jfe_success
      end if
   end function beta_value

   subroutine ols_fit(y, x, intercept, slope, residuals, status)
      real(dp), intent(in) :: y(:), x(:)
      real(dp), intent(out) :: intercept, slope
      real(dp), allocatable, intent(out) :: residuals(:)
      integer, intent(out) :: status
      real(dp), allocatable :: yc(:), xc(:)
      real(dp) :: mx, my, den

      call clean_pair(y, x, yc, xc)
      allocate(residuals(size(yc)))
      if (size(yc) < 2) then
         intercept = nan_dp()
         slope = nan_dp()
         residuals = nan_dp()
         status = jfe_insufficient_data
         return
      end if
      mx = sum(xc)/real(size(xc), dp)
      my = sum(yc)/real(size(yc), dp)
      den = sum((xc - mx)**2)
      if (den <= tiny(1.0_dp)) then
         intercept = nan_dp()
         slope = nan_dp()
         residuals = nan_dp()
         status = jfe_zero_denominator
         return
      end if
      slope = sum((xc - mx)*(yc - my))/den
      intercept = my - slope*mx
      residuals = yc - intercept - slope*xc
      status = jfe_success
   end subroutine ols_fit

   real(dp) function quantile_type1(x, probability, status) result(q)
      real(dp), intent(in) :: x(:), probability
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      real(dp) :: temp
      integer :: i, j, k

      call clean_vector(x, y)
      if (size(y) < 1 .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         q = nan_dp()
         if (present(status)) status = jfe_invalid_argument
         return
      end if
      do i = 2, size(y)
         temp = y(i)
         j = i - 1
         do while (j >= 1)
            if (y(j) <= temp) exit
            y(j + 1) = y(j)
            j = j - 1
         end do
         y(j + 1) = temp
      end do
      if (probability <= 0.0_dp) then
         k = 1
      else
         k = ceiling(probability*real(size(y), dp))
         k = max(1, min(k, size(y)))
      end if
      q = y(k)
      if (present(status)) status = jfe_success
   end function quantile_type1

   real(dp) function skewness_value(x, method, status) result(skew)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: method
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      real(dp) :: mu, m2, m3, s2
      integer :: n, meth

      meth = skew_moment
      if (present(method)) meth = method
      call clean_vector(x, y)
      n = size(y)
      if (n < 2) then
         skew = nan_dp()
         if (present(status)) status = jfe_insufficient_data
         return
      end if
      mu = sum(y)/real(n, dp)
      m2 = sum((y - mu)**2)/real(n, dp)
      m3 = sum((y - mu)**3)/real(n, dp)
      if (m2 <= tiny(1.0_dp)) then
         skew = nan_dp()
         if (present(status)) status = jfe_zero_denominator
         return
      end if
      select case (meth)
      case (skew_moment)
         skew = m3/m2**1.5_dp
      case (skew_fisher)
         if (n < 3) then
            skew = nan_dp()
            if (present(status)) status = jfe_insufficient_data
            return
         end if
         ! Mirrors the package's raw-moment Fisher definition.
         s2 = sum(y**2)/real(n, dp)
         if (s2 <= tiny(1.0_dp)) then
            skew = nan_dp()
            if (present(status)) status = jfe_zero_denominator
            return
         end if
         skew = sqrt(real(n*(n - 1), dp))/real(n - 2, dp) * &
            (sum(y**3)/real(n, dp))/s2**1.5_dp
      case (skew_sample)
         if (n < 3) then
            skew = nan_dp()
            if (present(status)) status = jfe_insufficient_data
            return
         end if
         skew = real(n*n, dp)/real((n - 1)*(n - 2), dp) * m3/m2**1.5_dp
      case default
         skew = nan_dp()
         if (present(status)) status = jfe_invalid_argument
         return
      end select
      if (present(status)) status = jfe_success
   end function skewness_value

   real(dp) function kurtosis_value(x, method, status) result(kurt)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: method
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      real(dp) :: mu, m2, m4, var_sample, raw2, raw4
      integer :: n, meth

      meth = kurt_excess
      if (present(method)) meth = method
      call clean_vector(x, y)
      n = size(y)
      if (n < 2) then
         kurt = nan_dp()
         if (present(status)) status = jfe_insufficient_data
         return
      end if
      mu = sum(y)/real(n, dp)
      m2 = sum((y - mu)**2)/real(n, dp)
      m4 = sum((y - mu)**4)/real(n, dp)
      if (m2 <= tiny(1.0_dp)) then
         kurt = nan_dp()
         if (present(status)) status = jfe_zero_denominator
         return
      end if
      select case (meth)
      case (kurt_excess)
         kurt = m4/m2**2 - 3.0_dp
      case (kurt_moment)
         kurt = m4/m2**2
      case (kurt_fisher)
         if (n < 4) then
            kurt = nan_dp()
            if (present(status)) status = jfe_insufficient_data
            return
         end if
         raw2 = sum(y**2)/real(n, dp)
         raw4 = sum(y**4)/real(n, dp)
         if (raw2 <= tiny(1.0_dp)) then
            kurt = nan_dp()
            if (present(status)) status = jfe_zero_denominator
            return
         end if
         kurt = real((n + 1)*(n - 1), dp) * &
            (raw4/raw2**2 - 3.0_dp*real(n - 1, dp)/real(n + 1, dp)) / &
            real((n - 2)*(n - 3), dp)
      case (kurt_sample, kurt_sample_excess)
         if (n < 4) then
            kurt = nan_dp()
            if (present(status)) status = jfe_insufficient_data
            return
         end if
         var_sample = sum((y - mu)**2)/real(n - 1, dp)
         kurt = sum((y - mu)**4)/var_sample**2 * &
            real(n*(n + 1), dp)/real((n - 1)*(n - 2)*(n - 3), dp)
         if (meth == kurt_sample_excess) then
            kurt = kurt - 3.0_dp*real((n - 1)*(n - 1), dp)/real((n - 2)*(n - 3), dp)
         end if
      case default
         kurt = nan_dp()
         if (present(status)) status = jfe_invalid_argument
         return
      end select
      if (present(status)) status = jfe_success
   end function kurtosis_value

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   real(dp) function annualized_return_value(x, scale, geometric, status) result(rann)
      real(dp), intent(in) :: x(:), scale
      logical, intent(in), optional :: geometric
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      logical :: use_geometric

      use_geometric = .true.
      if (present(geometric)) use_geometric = geometric
      call clean_vector(x, y)
      if (size(y) < 1 .or. scale <= 0.0_dp) then
         rann = nan_dp()
         if (present(status)) status = jfe_invalid_argument
         return
      end if
      if (use_geometric) then
         if (any(1.0_dp + y < 0.0_dp)) then
            rann = nan_dp()
            if (present(status)) status = jfe_invalid_argument
            return
         end if
         rann = product(1.0_dp + y)**(scale/real(size(y), dp)) - 1.0_dp
      else
         rann = sum(y)/real(size(y), dp)*scale
      end if
      if (present(status)) status = jfe_success
   end function annualized_return_value

   function drawdowns_value(x, geometric, status) result(dd)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: geometric
      integer, intent(out), optional :: status
      real(dp), allocatable :: dd(:)
      real(dp) :: cumulative, peak
      logical :: use_geometric
      integer :: i

      use_geometric = .true.
      if (present(geometric)) use_geometric = geometric
      allocate(dd(size(x)))
      cumulative = 1.0_dp
      peak = 1.0_dp
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) then
            dd(i) = nan_dp()
         else
            if (use_geometric) then
               cumulative = cumulative*(1.0_dp + x(i))
            else
               cumulative = cumulative + x(i)
            end if
            peak = max(peak, cumulative)
            if (abs(peak) <= tiny(1.0_dp)) then
               dd(i) = nan_dp()
            else
               dd(i) = cumulative/peak - 1.0_dp
            end if
         end if
      end do
      if (present(status)) status = jfe_success
   end function drawdowns_value

   function drawdown_peak_value(x, status) result(dd)
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: dd(:)
      real(dp) :: val
      integer :: i, j, peak

      allocate(dd(size(x)))
      dd = nan_dp()
      peak = 0
      do i = 1, size(x)
         if (.not. all(ieee_is_finite(x(peak + 1:i)))) cycle
         val = 1.0_dp
         do j = peak + 1, i
            val = val*(1.0_dp + x(j)/100.0_dp)
         end do
         if (val > 1.0_dp) then
            peak = i
            dd(i) = 0.0_dp
         else
            dd(i) = (val - 1.0_dp)*100.0_dp
         end if
      end do
      if (present(status)) status = jfe_success
   end function drawdown_peak_value

   real(dp) function max_drawdown_value(x, geometric, invert, status) result(mdd)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: geometric, invert
      integer, intent(out), optional :: status
      real(dp), allocatable :: dd(:), clean(:)
      logical :: inv

      inv = .true.
      if (present(invert)) inv = invert
      dd = drawdowns_value(x, geometric)
      call clean_vector(dd, clean)
      if (size(clean) < 1) then
         mdd = nan_dp()
         if (present(status)) status = jfe_insufficient_data
      else
         mdd = minval(clean)
         if (inv) mdd = -mdd
         if (present(status)) status = jfe_success
      end if
   end function max_drawdown_value

end module jfe_stats
