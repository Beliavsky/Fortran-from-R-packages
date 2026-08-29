! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_monitoring
   use r_kinds, only : dp
   use r_distributions, only : r_dnorm, r_pnorm
   use strucchange_tables, only : monitor_me_max, monitor_me_range
   use strucchange_tables, only : monitor_re_critvals
   use strucchange_utils, only : linear_interp
   implicit none
   private
   public :: log_plus
   public :: monitor_me_critical_value
   public :: monitor_ols_cusum_boundary
   public :: monitor_power_boundary
   public :: monitor_re_boundary
   public :: monitor_re_critical_value
   public :: mre_critical_value
   public :: pargmax_v
contains
   pure real(dp) function log_plus(x) result(value)
      real(dp), intent(in) :: x
      if (x <= exp(1.0_dp)) then
         value = 1.0_dp
      else
         value = log(x)
      end if
   end function log_plus

   real(dp) function mre_critical_value(alpha, tolerance) result(value)
      real(dp), intent(in) :: alpha
      real(dp), intent(in), optional :: tolerance
      real(dp) :: left, right, mid, fmid, tol
      integer :: iteration

      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = tolerance
      left = 0.0_dp
      right = 10.0_dp
      do iteration = 1, 200
         mid = 0.5_dp * (left + right)
         fmid = 2.0_dp * (r_pnorm(mid) - mid * r_dnorm(mid)) + alpha - 2.0_dp
         if (fmid < 0.0_dp) then
            left = mid
         else
            right = mid
         end if
         if (right - left <= tol * max(1.0_dp, mid)) exit
      end do
      value = 0.5_dp * (left + right)
   end function mre_critical_value

   real(dp) function monitor_me_critical_value(winsize, period, alpha, &
      functional, info) result(value)
      real(dp), intent(in) :: winsize, alpha
      integer, intent(in) :: period
      character(len = *), intent(in) :: functional
      integer, intent(out), optional :: info
      real(dp) :: pgrid(50), ygrid(50), target
      integer :: fidx, i, pidx, widx

      if (abs(winsize - 0.25_dp) < 1.0e-12_dp) then
         widx = 1
      else if (abs(winsize - 0.5_dp) < 1.0e-12_dp) then
         widx = 2
      else if (abs(winsize - 1.0_dp) < 1.0e-12_dp) then
         widx = 3
      else
         value = 0.0_dp
         if (present(info)) info = -1
         return
      end if
      select case (period)
      case (2)
         pidx = 1
      case (4)
         pidx = 2
      case (6)
         pidx = 3
      case (8)
         pidx = 4
      case (10)
         pidx = 5
      case default
         value = 0.0_dp
         if (present(info)) info = -2
         return
      end select
      select case (trim(functional))
      case ("max")
         fidx = 1
      case ("range")
         fidx = 2
      case default
         value = 0.0_dp
         if (present(info)) info = -3
         return
      end select
      do i = 1, 50
         pgrid(i) = 0.949_dp + 0.001_dp * real(i, dp)
         if (fidx == 1) then
            ygrid(i) = monitor_me_max(widx, pidx, i)
         else
            ygrid(i) = monitor_me_range(widx, pidx, i)
         end if
      end do
      target = 1.0_dp - alpha
      if (target < pgrid(1) .or. target > pgrid(50)) then
         value = 0.0_dp
         if (present(info)) info = -4
         return
      end if
      value = linear_interp(pgrid, ygrid, target)
      if (present(info)) info = 0
   end function monitor_me_critical_value

   real(dp) function monitor_re_critical_value(period, alpha, info) &
      result(value)
      integer, intent(in) :: period
      real(dp), intent(in) :: alpha
      integer, intent(out), optional :: info
      real(dp), parameter :: pgrid(2) = [0.9_dp, 0.95_dp]
      real(dp) :: target, ygrid(2)
      integer :: pidx

      select case (period)
      case (2)
         pidx = 1
      case (4)
         pidx = 2
      case (6)
         pidx = 3
      case (8)
         pidx = 4
      case (10)
         pidx = 5
      case default
         value = 0.0_dp
         if (present(info)) info = -1
         return
      end select
      target = 1.0_dp - alpha
      if (target < pgrid(1) .or. target > pgrid(2)) then
         value = 0.0_dp
         if (present(info)) info = -2
         return
      end if
      ygrid = monitor_re_critvals(pidx, :)
      value = linear_interp(pgrid, ygrid, target)
      if (present(info)) info = 0
   end function monitor_re_critical_value

   pure real(dp) function monitor_ols_cusum_boundary(k, history_size, &
      critical_value) result(value)
      integer, intent(in) :: k, history_size
      real(dp), intent(in) :: critical_value
      real(dp) :: x

      x = real(k, dp) / real(history_size, dp)
      if (x <= 1.0_dp) then
         value = 0.0_dp
      else
         value = sqrt(x * (x - 1.0_dp) * (critical_value ** 2 + &
            log(x / (x - 1.0_dp))))
      end if
   end function monitor_ols_cusum_boundary

   pure real(dp) function monitor_re_boundary(k, history_size, &
      critical_value) result(value)
      integer, intent(in) :: k, history_size
      real(dp), intent(in) :: critical_value
      real(dp) :: x

      x = real(k, dp) / real(history_size, dp)
      if (x <= 1.0_dp) then
         value = 0.0_dp
      else
         value = sqrt(x * (x - 1.0_dp) * (critical_value ** 2 + &
            log(x / (x - 1.0_dp))))
      end if
   end function monitor_re_boundary

   pure real(dp) function monitor_power_boundary(k, history_size, &
      critical_value) result(value)
      integer, intent(in) :: k, history_size
      real(dp), intent(in) :: critical_value

      value = critical_value * sqrt(2.0_dp * &
         log_plus(real(k, dp) / real(history_size, dp)))
   end function monitor_power_boundary

   real(dp) function pargmax_v(x, xi, phi1, phi2) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: xi, phi1, phi2
      real(dp) :: ax, frac, phi, xi_use, phi1_use, phi2_use
      real(dp), parameter :: two_pi = 6.2831853071795864769_dp

      xi_use = 1.0_dp
      phi1_use = 1.0_dp
      phi2_use = 1.0_dp
      if (present(xi)) xi_use = xi
      if (present(phi1)) phi1_use = phi1
      if (present(phi2)) phi2_use = phi2
      phi = xi_use * (phi2_use / phi1_use) ** 2
      ax = abs(x)
      if (x < 0.0_dp) then
         frac = xi_use / phi
         value = -exp(0.5_dp * log(ax) - ax / 8.0_dp - 0.5_dp * log(two_pi)) - &
            (phi / xi_use * (phi + 2.0_dp * xi_use) / (phi + xi_use)) * &
            exp(frac * (1.0_dp + frac) * ax / 2.0_dp + &
            r_pnorm(-(0.5_dp + frac) * sqrt(ax), log_probability = .true.)) + &
            exp(log(ax / 2.0_dp - 2.0_dp + &
            (phi + 2.0_dp * xi_use) ** 2 / ((phi + xi_use) * xi_use)) + &
            r_pnorm(-sqrt(ax) / 2.0_dp, log_probability = .true.))
      else
         frac = xi_use ** 2 / phi
         value = 1.0_dp + sqrt(frac) * &
            exp(0.5_dp * log(max(ax, tiny(1.0_dp))) - frac * ax / 8.0_dp - &
            0.5_dp * log(two_pi)) + &
            (xi_use / phi * (2.0_dp * phi + xi_use) / (phi + xi_use)) * &
            exp((phi + xi_use) * ax / 2.0_dp + &
            r_pnorm(-(phi + xi_use / 2.0_dp) / sqrt(phi) * sqrt(ax), &
            log_probability = .true.)) - &
            exp(log((2.0_dp * phi + xi_use) ** 2 / ((phi + xi_use) * phi) - &
            2.0_dp + frac * ax / 2.0_dp) + &
            r_pnorm(-sqrt(frac) * sqrt(ax) / 2.0_dp, log_probability = .true.))
      end if
   end function pargmax_v
end module strucchange_monitoring
