! Numerical helpers for argus-fortran.
! SPDX-License-Identifier: GPL-2.0-or-later
module argus_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_negative_inf
   use argus_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: a32 = 1.5_dp
   real(dp), parameter :: log_gamma_3half = 0.5_dp*log(acos(-1.0_dp)) - log(2.0_dp)
   real(dp), parameter :: log_gamma_5half = log(3.0_dp) + 0.5_dp*log(acos(-1.0_dp)) - log(4.0_dp)
   real(dp), parameter :: two_over_sqrt_pi = 2.0_dp/sqrt(acos(-1.0_dp))

   public :: gamma_p_3half, log_gamma_p_3half
   public :: inv_gamma_p_3half_log
   public :: log1p_stable, expm1_stable, log1mexp
   public :: nan_dp, neg_inf_dp

contains

   pure real(dp) function nan_dp() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure real(dp) function neg_inf_dp() result(x)
      x = ieee_value(0.0_dp, ieee_negative_inf)
   end function neg_inf_dp

   pure real(dp) function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: term, sum
      integer :: k

      if (x <= -1.0_dp) then
         if (x < -1.0_dp) then
            y = nan_dp()
         else
            y = neg_inf_dp()
         end if
      else if (abs(x) > 1.0e-4_dp) then
         y = log(1.0_dp + x)
      else
         term = x
         sum = term
         do k = 2, 30
            term = term*x
            if (mod(k,2) == 0) then
               sum = sum - term/real(k,dp)
            else
               sum = sum + term/real(k,dp)
            end if
            if (abs(term/real(k,dp)) < epsilon(1.0_dp)*max(1.0_dp,abs(sum))) exit
         end do
         y = sum
      end if
   end function log1p_stable

   pure real(dp) function expm1_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: term, sum
      integer :: k

      if (abs(x) > 1.0e-5_dp) then
         y = exp(x) - 1.0_dp
      else
         term = x
         sum = term
         do k = 2, 30
            term = term*x/real(k,dp)
            sum = sum + term
            if (abs(term) < epsilon(1.0_dp)*max(1.0_dp,abs(sum))) exit
         end do
         y = sum
      end if
   end function expm1_stable

   pure real(dp) function log1mexp(logx) result(y)
      real(dp), intent(in) :: logx
      real(dp), parameter :: log_half = -0.693147180559945309417232121458176568_dp

      if (logx >= 0.0_dp) then
         if (logx > 0.0_dp) then
            y = nan_dp()
         else
            y = neg_inf_dp()
         end if
      else if (logx < log_half) then
         y = log1p_stable(-exp(logx))
      else
         y = log(-expm1_stable(logx))
      end if
   end function log1mexp

   pure real(dp) function log_gamma_p_3half(x) result(logp)
      real(dp), intent(in) :: x
      real(dp) :: ap, del, s, q, rootx
      integer :: n

      if (x <= 0.0_dp) then
         if (x < 0.0_dp) then
            logp = nan_dp()
         else
            logp = neg_inf_dp()
         end if
         return
      end if

      ! For small x, the usual lower-incomplete-gamma series avoids
      ! catastrophic cancellation in 1-Q(a,x).
      if (x < 1.0_dp) then
         ap = a32
         s = 1.0_dp/a32
         del = s
         do n = 1, 10000
            ap = ap + 1.0_dp
            del = del*x/ap
            s = s + del
            if (abs(del) <= epsilon(1.0_dp)*abs(s)) exit
         end do
         logp = -x + a32*log(x) - log_gamma_3half + log(s)
      else
         ! Q(3/2,x) = erfc(sqrt(x)) + 2 sqrt(x) exp(-x)/sqrt(pi).
         rootx = sqrt(x)
         q = erfc(rootx) + two_over_sqrt_pi*rootx*exp(-x)
         if (q <= 0.0_dp) then
            logp = 0.0_dp
         else if (q >= 1.0_dp) then
            ! This branch is only a roundoff guard; x>=1 implies Q<1.
            logp = neg_inf_dp()
         else
            logp = log1p_stable(-q)
         end if
      end if
   end function log_gamma_p_3half

   pure real(dp) function gamma_p_3half(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: lp

      lp = log_gamma_p_3half(x)
      p = exp(lp)
   end function gamma_p_3half

   pure real(dp) function inv_gamma_p_3half_log(log_target, upper_bound) result(x)
      real(dp), intent(in) :: log_target, upper_bound
      real(dp) :: lo, hi, y, ynew, lp, lpdf, deriv, approx
      real(dp) :: lp_upper, scale
      integer :: iter

      if (upper_bound < 0.0_dp .or. log_target > 0.0_dp) then
         x = nan_dp()
         return
      end if
      if (log_target < -huge(1.0_dp)) then
         x = 0.0_dp
         return
      end if

      lp_upper = log_gamma_p_3half(upper_bound)
      if (log_target >= lp_upper - 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lp_upper))) then
         x = upper_bound
         return
      end if

      ! Small-P asymptotic inversion: P(3/2,x) ~ x^(3/2)/Gamma(5/2).
      approx = exp((log_target + log_gamma_5half)/a32)
      if (.not. (approx > 0.0_dp)) approx = min(upper_bound, tiny(1.0_dp)**(2.0_dp/3.0_dp))

      lo = 0.0_dp
      hi = min(upper_bound, max(1.0_dp, 2.0_dp*approx + 0.25_dp))
      do while (hi < upper_bound)
         if (log_gamma_p_3half(hi) >= log_target) exit
         hi = min(upper_bound, 2.0_dp*hi + 1.0_dp)
      end do

      y = min(hi, max(approx, 0.5_dp*hi))
      do iter = 1, 100
         lp = log_gamma_p_3half(y)
         if (lp < log_target) then
            lo = y
         else
            hi = y
         end if

         scale = max(1.0_dp, abs(y))
         if (hi-lo <= 8.0_dp*epsilon(1.0_dp)*scale) exit

         if (y > 0.0_dp .and. lp > neg_inf_dp()) then
            lpdf = 0.5_dp*log(y) - y - log_gamma_3half
            if (lpdf-lp < log(huge(1.0_dp))) then
               deriv = exp(lpdf-lp)
            else
               deriv = huge(1.0_dp)
            end if
            if (deriv > 0.0_dp) then
               ynew = y - (lp-log_target)/deriv
            else
               ynew = 0.5_dp*(lo+hi)
            end if
         else
            ynew = 0.5_dp*(lo+hi)
         end if

         if (.not. (ynew > lo .and. ynew < hi)) ynew = 0.5_dp*(lo+hi)
         y = ynew
      end do
      x = 0.5_dp*(lo+hi)
   end function inv_gamma_p_3half_log

end module argus_special
