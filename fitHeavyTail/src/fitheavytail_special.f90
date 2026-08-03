! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_special
   use fitheavytail_kinds, only: dp
   implicit none
   private
   public :: digamma_dp, log_bessel_k, bessel_k_ratio, bessel_order_derivative

contains

   elemental function digamma_dp(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, y, inv, inv2

      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      y = x
      value = 0.0_dp
      do while (y < 8.0_dp)
         value = value - 1.0_dp/y
         y = y + 1.0_dp
      end do
      inv = 1.0_dp/y
      inv2 = inv*inv
      value = value + log(y) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - &
         inv2*(1.0_dp/120.0_dp - inv2*(1.0_dp/252.0_dp - &
         inv2*(1.0_dp/240.0_dp - inv2*(5.0_dp/660.0_dp)))))
   end function digamma_dp

   function log_bessel_k(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value
      real(dp) :: order, frac, log_prev, log_curr, log_next, coefficient
      integer :: m, k

      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      order = abs(nu)
      m = floor(order)
      frac = order - real(m,dp)
      log_prev = log_bessel_k_base(x, frac)
      if (m == 0) then
         value = log_prev
         return
      end if
      log_curr = log_bessel_k_base(x, frac + 1.0_dp)
      if (m == 1) then
         value = log_curr
         return
      end if
      do k = 1, m-1
         coefficient = 2.0_dp*(frac + real(k,dp))/x
         log_next = log_add_exp(log_prev, log(coefficient)+log_curr)
         log_prev = log_curr
         log_curr = log_next
      end do
      value = log_curr
   end function log_bessel_k

   function bessel_k_ratio(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value, d
      d = log_bessel_k(x, nu+1.0_dp) - log_bessel_k(x, nu)
      if (d > log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else if (d < log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(d)
      end if
   end function bessel_k_ratio

   function bessel_order_derivative(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value, h
      h = 2.0e-4_dp*(1.0_dp + abs(nu))
      value = (log_bessel_k(x,nu+h)-log_bessel_k(x,nu-h))/(2.0_dp*h)
   end function bessel_order_derivative

   function log_bessel_k_base(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value, upper, integral, fa, fb, fm, whole
      integer :: depth

      upper = max(20.0_dp, acosh(1.0_dp + 60.0_dp/max(x,1.0e-300_dp)))
      upper = min(40.0_dp, upper)
      fa = scaled_integrand(0.0_dp, x, nu)
      fb = scaled_integrand(upper, x, nu)
      fm = scaled_integrand(0.5_dp*upper, x, nu)
      whole = upper*(fa + 4.0_dp*fm + fb)/6.0_dp
      depth = 22
      integral = adaptive_simpson(0.0_dp, upper, fa, fm, fb, whole, &
         1.0e-11_dp*max(1.0_dp,abs(whole)), depth, x, nu)
      value = -x + log(max(integral,tiny(1.0_dp)))
   end function log_bessel_k_base

   recursive function adaptive_simpson(a, b, fa, fm, fb, whole, tol, depth, x, nu) result(value)
      real(dp), intent(in) :: a, b, fa, fm, fb, whole, tol, x, nu
      integer, intent(in) :: depth
      real(dp) :: value, left_mid, right_mid, flm, frm, left, right, mid

      mid = 0.5_dp*(a+b)
      left_mid = 0.5_dp*(a+mid)
      right_mid = 0.5_dp*(mid+b)
      flm = scaled_integrand(left_mid,x,nu)
      frm = scaled_integrand(right_mid,x,nu)
      left = (mid-a)*(fa + 4.0_dp*flm + fm)/6.0_dp
      right = (b-mid)*(fm + 4.0_dp*frm + fb)/6.0_dp
      if (depth <= 0 .or. abs(left+right-whole) <= 15.0_dp*tol) then
         value = left + right + (left+right-whole)/15.0_dp
      else
         value = adaptive_simpson(a,mid,fa,flm,fm,left,0.5_dp*tol,depth-1,x,nu) + &
                 adaptive_simpson(mid,b,fm,frm,fb,right,0.5_dp*tol,depth-1,x,nu)
      end if
   end function adaptive_simpson

   function scaled_integrand(t, x, nu) result(value)
      real(dp), intent(in) :: t, x, nu
      real(dp) :: value, exponent
      exponent = -x*(cosh(t)-1.0_dp) + log_cosh(nu*t)
      if (exponent < log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(exponent)
      end if
   end function scaled_integrand

   elemental function log_cosh(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, ax
      ax = abs(x)
      if (ax < 20.0_dp) then
         value = log(cosh(x))
      else
         value = ax - log(2.0_dp) + log(1.0_dp + exp(-2.0_dp*ax))
      end if
   end function log_cosh

   elemental function log_add_exp(a,b) result(value)
      real(dp), intent(in) :: a,b
      real(dp) :: value, m
      m = max(a,b)
      value = m + log(exp(a-m)+exp(b-m))
   end function log_add_exp

end module fitheavytail_special
