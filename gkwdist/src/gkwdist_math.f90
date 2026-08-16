! SPDX-License-Identifier: MIT
module gkwdist_math
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, ieee_is_finite
   use gkwdist_kinds, only : dp
   implicit none
   private
   public :: log1mexp, safe_log, safe_exp, safe_pow, log_beta, digamma_fn, trigamma_fn
   public :: beta_cdf, beta_quantile, beta_logpdf, nan_dp, posinf_dp, neginf_dp, finite_dp
   public :: log1p_stable, expm1_stable

contains

   pure function expm1_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (abs(x) < 1.0e-5_dp) then
         y = x*(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x*(1.0_dp/24.0_dp + x/120.0_dp))))
      else
         y = exp(x)-1.0_dp
      end if
   end function expm1_stable

   pure function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y, z
      if (x <= -1.0_dp) then
         if (x == -1.0_dp) then
            y = neginf_dp()
         else
            y = nan_dp()
         end if
      else if (abs(x) < 1.0e-4_dp) then
         z=x
         y=z-z*z/2.0_dp+z**3/3.0_dp-z**4/4.0_dp+z**5/5.0_dp-z**6/6.0_dp
      else
         y=log(1.0_dp+x)
      end if
   end function log1p_stable

   pure function nan_dp() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure function posinf_dp() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function posinf_dp

   pure function neginf_dp() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_negative_inf)
   end function neginf_dp

   pure logical function finite_dp(x)
      real(dp), intent(in) :: x
      finite_dp = ieee_is_finite(x)
   end function finite_dp

   pure function log1mexp(u) result(y)
      real(dp), intent(in) :: u
      real(dp) :: y
      real(dp), parameter :: crossover = -0.6931471805599453094172321214581766_dp
      if (u > 0.0_dp) then
         y = nan_dp()
      else if (u == 0.0_dp) then
         y = neginf_dp()
      else if (u > -1.0e-14_dp) then
         y = log(-u) - 0.5_dp*u
      else if (u > crossover) then
         y = log(-expm1_stable(u))
      else
         y = log1p_stable(-exp(u))
      end if
   end function log1mexp

   pure function safe_log(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (x < 0.0_dp) then
         y = nan_dp()
      else if (x == 0.0_dp) then
         y = neginf_dp()
      else
         y = log(x)
      end if
   end function safe_log

   pure function safe_exp(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      real(dp), parameter :: log_huge = log(huge(1.0_dp))
      real(dp), parameter :: log_tiny = log(tiny(1.0_dp))
      if (x > log_huge) then
         y = posinf_dp()
      else if (x < log_tiny - 10.0_dp) then
         y = 0.0_dp
      else
         y = exp(x)
      end if
   end function safe_exp

   pure function safe_pow(x, y) result(z)
      real(dp), intent(in) :: x, y
      real(dp) :: z, yi
      if (.not. finite_dp(x) .or. .not. finite_dp(y)) then
         z = nan_dp()
      else if (x == 0.0_dp) then
         if (y > 0.0_dp) then
            z = 0.0_dp
         else if (y == 0.0_dp) then
            z = 1.0_dp
         else
            z = posinf_dp()
         end if
      else if (x < 0.0_dp) then
         yi = anint(y)
         if (abs(y - yi) > 1.0e-12_dp) then
            z = nan_dp()
         else
            z = x**int(yi)
         end if
      else
         z = safe_exp(y*log(x))
      end if
   end function safe_pow

   pure function log_beta(a, b) result(v)
      real(dp), intent(in) :: a, b
      real(dp) :: v
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = nan_dp()
      else
         v = log_gamma(a) + log_gamma(b) - log_gamma(a+b)
      end if
   end function log_beta

   pure function digamma_fn(xin) result(v)
      real(dp), intent(in) :: xin
      real(dp) :: v, x, inv, inv2
      if (xin <= 0.0_dp .or. .not. finite_dp(xin)) then
         v = nan_dp()
         return
      end if
      x = xin
      v = 0.0_dp
      do while (x < 8.0_dp)
         v = v - 1.0_dp/x
         x = x + 1.0_dp
      end do
      inv = 1.0_dp/x
      inv2 = inv*inv
      v = v + log(x) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
         inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp - inv2*(5.0_dp/660.0_dp)))))
   end function digamma_fn

   pure function trigamma_fn(xin) result(v)
      real(dp), intent(in) :: xin
      real(dp) :: v, x, inv, inv2
      if (xin <= 0.0_dp .or. .not. finite_dp(xin)) then
         v = nan_dp()
         return
      end if
      x = xin
      v = 0.0_dp
      do while (x < 8.0_dp)
         v = v + 1.0_dp/(x*x)
         x = x + 1.0_dp
      end do
      inv = 1.0_dp/x
      inv2 = inv*inv
      v = v + inv + 0.5_dp*inv2 + inv*inv2/6.0_dp - inv*inv2*inv2/30.0_dp + &
         inv*inv2*inv2*inv2/42.0_dp - inv*inv2*inv2*inv2*inv2/30.0_dp + &
         5.0_dp*inv*inv2*inv2*inv2*inv2*inv2/66.0_dp
   end function trigamma_fn

   pure function betacf(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      real(dp) :: cf
      integer, parameter :: maxit = 300
      real(dp), parameter :: eps = 3.0e-15_dp, fpmin = tiny(1.0_dp)/eps
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap
      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      cf = h
   end function betacf

   pure function beta_cdf(x, a, b, lower_tail) result(p)
      real(dp), intent(in) :: x, a, b
      logical, intent(in), optional :: lower_tail
      real(dp) :: p, bt, lower, upper
      logical :: lower_flag
      lower_flag = .true.
      if (present(lower_tail)) lower_flag = lower_tail
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         p = nan_dp()
         return
      end if
      if (x <= 0.0_dp) then
         p = merge(0.0_dp, 1.0_dp, lower_flag)
         return
      else if (x >= 1.0_dp) then
         p = merge(1.0_dp, 0.0_dp, lower_flag)
         return
      end if
      bt = exp(a*log(x) + b*log1p_stable(-x) - log_beta(a,b))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         lower = bt*betacf(a,b,x)/a
         lower = max(0.0_dp, min(1.0_dp, lower))
         upper = 1.0_dp - lower
      else
         upper = bt*betacf(b,a,1.0_dp-x)/b
         upper = max(0.0_dp, min(1.0_dp, upper))
         lower = 1.0_dp - upper
      end if
      p = merge(lower, upper, lower_flag)
   end function beta_cdf

   pure function beta_logpdf(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      real(dp) :: v
      if (x <= 0.0_dp .or. x >= 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         v = neginf_dp()
      else
         v = (a-1.0_dp)*log(x) + (b-1.0_dp)*log1p_stable(-x) - log_beta(a,b)
      end if
   end function beta_logpdf

   pure function beta_quantile(p, a, b, lower_tail, log_p) result(x)
      real(dp), intent(in) :: p, a, b
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: x, pp, lo, hi, fx, dens, xn
      logical :: lower_flag, log_flag
      integer :: it
      lower_flag = .true.; log_flag = .false.
      if (present(lower_tail)) lower_flag = lower_tail
      if (present(log_p)) log_flag = log_p
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         x = nan_dp(); return
      end if
      if (log_flag) then
         if (p > 0.0_dp) then
            x = nan_dp(); return
         end if
         pp = exp(p)
      else
         pp = p
      end if
      if (.not. lower_flag) pp = 1.0_dp - pp
      if (pp <= 0.0_dp) then
         x = 0.0_dp; return
      else if (pp >= 1.0_dp) then
         x = 1.0_dp; return
      end if
      lo = 0.0_dp; hi = 1.0_dp
      x = a/(a+b)
      x = max(1.0e-12_dp, min(1.0_dp-1.0e-12_dp, x))
      do it = 1, 100
         fx = beta_cdf(x,a,b,.true.) - pp
         if (fx > 0.0_dp) then
            hi = x
         else
            lo = x
         end if
         dens = exp(beta_logpdf(x,a,b))
         if (dens > 0.0_dp .and. finite_dp(dens)) then
            xn = x - fx/dens
         else
            xn = 0.5_dp*(lo+hi)
         end if
         if (xn <= lo .or. xn >= hi .or. .not. finite_dp(xn)) xn = 0.5_dp*(lo+hi)
         if (abs(xn-x) <= 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x))) then
            x = xn
            exit
         end if
         x = xn
      end do
      x = max(0.0_dp, min(1.0_dp, x))
   end function beta_quantile

end module gkwdist_math
