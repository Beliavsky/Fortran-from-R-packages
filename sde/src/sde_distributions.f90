! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_distributions
   use sde_kinds, only : dp
   use sde_special, only : chi_square_cdf, nan_dp, log_sum_exp
   implicit none
   private

   public :: noncentral_chi_square_pdf
   public :: noncentral_chi_square_logpdf
   public :: noncentral_chi_square_cdf
   public :: noncentral_chi_square_quantile

contains

   pure function central_chi_square_logpdf(x, df) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: df
      real(dp) :: value

      if (df <= 0.0_dp) then
         value = nan_dp()
      else if (x < 0.0_dp) then
         value = -huge(1.0_dp)
      else if (x <= 0.0_dp) then
         if (df < 2.0_dp) then
            value = huge(1.0_dp)
         else if (df <= 2.0_dp) then
            value = -log(2.0_dp)
         else
            value = -huge(1.0_dp)
         end if
      else
         value = (0.5_dp*df-1.0_dp)*log(x)-0.5_dp*x-0.5_dp*df*log(2.0_dp)-log_gamma(0.5_dp*df)
      end if
   end function central_chi_square_logpdf

   pure subroutine poisson_window(ncp, k_low, k_high)
      real(dp), intent(in) :: ncp
      integer, intent(out) :: k_low
      integer, intent(out) :: k_high
      real(dp) :: mean_k, spread

      mean_k = 0.5_dp*ncp
      spread = 12.0_dp*sqrt(mean_k+1.0_dp)+20.0_dp
      k_low = max(0, floor(mean_k-spread))
      k_high = max(k_low+20, ceiling(mean_k+spread))
   end subroutine poisson_window

   pure function noncentral_chi_square_logpdf(x, df, ncp) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: df
      real(dp), intent(in) :: ncp
      real(dp) :: value
      real(dp), allocatable :: terms(:)
      real(dp) :: lambda
      integer :: k, k_low, k_high, index

      if (df <= 0.0_dp .or. ncp < 0.0_dp) then
         value = nan_dp()
         return
      end if
      if (ncp <= 0.0_dp) then
         value = central_chi_square_logpdf(x, df)
         return
      end if
      if (x < 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      lambda = 0.5_dp*ncp
      call poisson_window(ncp, k_low, k_high)
      allocate(terms(k_high-k_low+1))
      do k = k_low, k_high
         index = k-k_low+1
         terms(index) = -lambda+real(k, dp)*log(lambda)-log_gamma(real(k+1, dp))+ &
            central_chi_square_logpdf(x, df+2.0_dp*real(k, dp))
      end do
      value = log_sum_exp(terms)
   end function noncentral_chi_square_logpdf

   pure function noncentral_chi_square_pdf(x, df, ncp) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: df
      real(dp), intent(in) :: ncp
      real(dp) :: value
      real(dp) :: log_value

      log_value = noncentral_chi_square_logpdf(x, df, ncp)
      if (log_value >= log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else if (log_value <= log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(log_value)
      end if
   end function noncentral_chi_square_pdf

   pure function noncentral_chi_square_cdf(x, df, ncp, lower_tail) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: df
      real(dp), intent(in) :: ncp
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      real(dp) :: lambda, maximum_log_weight, weight, sum_weight, sum_value
      real(dp), allocatable :: log_weights(:)
      integer :: k, k_low, k_high, index
      logical :: lower

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (df <= 0.0_dp .or. ncp < 0.0_dp) then
         value = nan_dp()
         return
      end if
      if (x <= 0.0_dp) then
         if (lower) then
            value = 0.0_dp
         else
            value = 1.0_dp
         end if
         return
      end if
      if (ncp <= 0.0_dp) then
         value = chi_square_cdf(x, df, lower)
         return
      end if

      lambda = 0.5_dp*ncp
      call poisson_window(ncp, k_low, k_high)
      allocate(log_weights(k_high-k_low+1))
      do k = k_low, k_high
         index = k-k_low+1
         log_weights(index) = -lambda+real(k, dp)*log(lambda)-log_gamma(real(k+1, dp))
      end do
      maximum_log_weight = maxval(log_weights)
      sum_weight = 0.0_dp
      sum_value = 0.0_dp
      do k = k_low, k_high
         index = k-k_low+1
         weight = exp(log_weights(index)-maximum_log_weight)
         sum_weight = sum_weight+weight
         sum_value = sum_value+weight*chi_square_cdf(x, df+2.0_dp*real(k, dp))
      end do
      value = sum_value/sum_weight
      if (.not. lower) value = 1.0_dp-value
      value = max(0.0_dp, min(1.0_dp, value))
   end function noncentral_chi_square_cdf

   pure function noncentral_chi_square_quantile(p, df, ncp, lower_tail) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in) :: df
      real(dp), intent(in) :: ncp
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      real(dp) :: prob, lower_x, upper_x, midpoint, mean_value, variance_value
      integer :: iter
      logical :: lower

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (df <= 0.0_dp .or. ncp < 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         value = nan_dp()
         return
      end if
      prob = p
      if (.not. lower) prob = 1.0_dp-p
      if (prob <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (prob >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if

      mean_value = df+ncp
      variance_value = 2.0_dp*(df+2.0_dp*ncp)
      lower_x = 0.0_dp
      upper_x = max(1.0_dp, mean_value+12.0_dp*sqrt(variance_value)+20.0_dp)
      do while (noncentral_chi_square_cdf(upper_x, df, ncp) < prob)
         upper_x = 2.0_dp*upper_x
      end do
      do iter = 1, 220
         midpoint = 0.5_dp*(lower_x+upper_x)
         if (noncentral_chi_square_cdf(midpoint, df, ncp) < prob) then
            lower_x = midpoint
         else
            upper_x = midpoint
         end if
         if (upper_x-lower_x <= 2.0e-12_dp*max(1.0_dp, midpoint)) exit
      end do
      value = 0.5_dp*(lower_x+upper_x)
   end function noncentral_chi_square_quantile

end module sde_distributions
