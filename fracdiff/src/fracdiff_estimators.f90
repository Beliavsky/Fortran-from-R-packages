! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_estimators
   use fracdiff_kinds, only : dp, pi_dp
   use fracdiff_status, only : fd_ok, fd_invalid_input, fd_insufficient_data, &
                               fracdiff_status_message
   use fracdiff_types, only : fractional_d_estimate
   use r_descriptive, only : r_mean
   implicit none
   private

   public :: fd_gph, fd_sperio

contains

   function fd_gph(x, bandwidth_exponent) result(estimate)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: bandwidth_exponent
      type(fractional_d_estimate) :: estimate
      real(dp) :: exponent

      exponent = 0.5_dp
      if (present(bandwidth_exponent)) exponent = bandwidth_exponent
      call estimate_from_periodogram(x, exponent, .false., 0.9_dp, estimate)
   end function fd_gph

   function fd_sperio(x, bandwidth_exponent, beta) result(estimate)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: bandwidth_exponent, beta
      type(fractional_d_estimate) :: estimate
      real(dp) :: exponent, beta_value

      exponent = 0.5_dp
      beta_value = 0.9_dp
      if (present(bandwidth_exponent)) exponent = bandwidth_exponent
      if (present(beta)) beta_value = beta
      call estimate_from_periodogram(x, exponent, .true., beta_value, estimate)
   end function fd_sperio

   subroutine estimate_from_periodogram(x, bandwidth_exponent, tapered, beta, estimate)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: bandwidth_exponent, beta
      logical, intent(in) :: tapered
      type(fractional_d_estimate), intent(out) :: estimate

      real(dp), allocatable :: centered(:), autocovariance(:), weights(:)
      real(dp), allocatable :: xreg(:), yreg(:)
      real(dp) :: mean_x, variance_x, frequency, periodogram, a_k
      real(dp) :: xbar, ybar, sxx, sxy, residual_sum, slope, variance_d
      integer :: n, g, m, m2, k, i, n_positive

      n = size(x)
      estimate%status = fd_ok
      if (n < 8 .or. bandwidth_exponent <= 0.0_dp .or. bandwidth_exponent >= 1.0_dp) then
         estimate%status = fd_invalid_input
         estimate%message = fracdiff_status_message(estimate%status)
         return
      end if
      if (tapered .and. (beta <= 0.0_dp .or. beta > 1.0_dp)) then
         estimate%status = fd_invalid_input
         estimate%message = fracdiff_status_message(estimate%status)
         return
      end if

      g = int(real(n,dp)**bandwidth_exponent)
      if (g < 2) then
         estimate%status = fd_insufficient_data
         estimate%message = fracdiff_status_message(estimate%status)
         return
      end if

      allocate(centered(n), autocovariance(n-1), weights(n-1), xreg(g), yreg(g))
      mean_x = r_mean(x)
      centered = x - mean_x
      variance_x = dot_product(centered,centered)/real(n,dp)
      do k = 1, n - 1
         autocovariance(k) = dot_product(centered(1:n-k), centered(k+1:n))/real(n,dp)
      end do

      weights = 1.0_dp
      m = 0
      if (tapered) then
         m = int(real(n,dp)**beta)
         m = max(1,m)
         m2 = m/2
         do k = 1, n - 1
            a_k = real(k,dp)/real(m,dp)
            if (k <= m2) then
               weights(k) = 1.0_dp - 6.0_dp*a_k*a_k*(1.0_dp - a_k)
            else if (k <= m) then
               weights(k) = 2.0_dp*(1.0_dp - a_k)**3
            else
               weights(k) = 0.0_dp
            end if
         end do
      end if

      n_positive = 0
      do i = 1, g
         frequency = 2.0_dp*pi_dp*real(i,dp)/real(n,dp)
         periodogram = variance_x
         do k = 1, n - 1
            periodogram = periodogram + 2.0_dp*autocovariance(k)*weights(k)* &
               cos(frequency*real(k,dp))
         end do
         if (periodogram > 0.0_dp) then
            n_positive = n_positive + 1
            yreg(n_positive) = log(periodogram/(2.0_dp*pi_dp))
            xreg(n_positive) = 2.0_dp*log(2.0_dp*sin(frequency/2.0_dp))
         end if
      end do

      if (n_positive < 2) then
         estimate%status = fd_insufficient_data
         estimate%message = "fewer than two positive periodogram ordinates"
         return
      end if
      xbar = r_mean(xreg(1:n_positive))
      ybar = r_mean(yreg(1:n_positive))
      sxx = sum((xreg(1:n_positive)-xbar)**2)
      sxy = dot_product(xreg(1:n_positive)-xbar, yreg(1:n_positive)-ybar)
      if (sxx <= 0.0_dp) then
         estimate%status = fd_insufficient_data
         estimate%message = fracdiff_status_message(estimate%status)
         return
      end if
      slope = sxy/sxx
      residual_sum = sum((yreg(1:n_positive) - (ybar + slope*(xreg(1:n_positive)-xbar)))**2)

      estimate%d = -slope
      if (tapered) then
         variance_d = (0.539285_dp*real(m,dp)/real(n,dp))/sxx
      else
         variance_d = pi_dp*pi_dp/(6.0_dp*sxx)
      end if
      estimate%sd_asymptotic = sqrt(max(variance_d,0.0_dp))
      estimate%sd_regression = sqrt(max(residual_sum/(real(g-1,dp)*sxx),0.0_dp))
      estimate%message = "ok"
   end subroutine estimate_from_periodogram

end module fracdiff_estimators
