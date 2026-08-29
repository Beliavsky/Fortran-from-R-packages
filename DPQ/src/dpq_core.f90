! Core numerical helpers translated from DPQ.
! SPDX-License-Identifier: GPL-2.0-or-later
module dpq_core
   use r_compat, only: dp, r_lgamma, r_gamma, r_psigamma
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
      ieee_positive_inf, ieee_negative_inf, ieee_is_finite
   implicit none
   private

   real(dp), parameter, public :: m_ln2 = 0.693147180559945309417232121458176568_dp
   real(dp), parameter, public :: m_sqrt2 = 1.41421356237309504880168872420969808_dp
   real(dp), parameter, public :: m_ln_sqrt_2pi = 0.91893853320467274178032973640561764_dp
   real(dp), parameter, public :: g_half = 1.77245385090551602729816748334114518_dp
   real(dp), parameter, public :: m_cutoff = m_ln2 * real(digits(1.0_dp), dp)
   real(dp), parameter, public :: m_minexp = log(tiny(1.0_dp))

   public :: d_0, d_1, dt_0, dt_1, d_lval, d_cval, d_val, d_qiv, d_exp, d_log, d_clog
   public :: d_lexp, dt_val, dt_cval, dt_qiv, dt_civ, dt_log, dt_clog, dt_log_known
   public :: log1p_dp, expm1_dp, log1mexp, log1pexp, logspace_add, logspace_sub
   public :: pow1p, pow_di, dpq_pow, expm1x, expm1x_tser, rexpm1, rlog1, log4p1p, p1l1, p1l1p, p1l1ser
   public :: logcf, log1pmx, lgamma1p, lgamma1p_series, logr
   public :: lsum, lssum, chebyshev_eval, chebyshev_poly, chebyshev_nc
   public :: dpq_frexp, dpq_ldexp, dpq_modf, dpsifn_scalar
   public :: clamp01, prob_from_input, prob_output, lower_prob_log

contains

   pure elemental real(dp) function log1p_dp(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (x == -1.0_dp) then
         v = ieee_value(0.0_dp, ieee_negative_inf)
      else if (x < -1.0_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (abs(x) > 0.5_dp) then
         v = log(1.0_dp + x)
      else
         y = 1.0_dp + x
         if (y == 1.0_dp) then
            v = x
         else
            v = log(y) - ((y - 1.0_dp) - x)/y
         end if
      end if
   end function log1p_dp

   pure elemental real(dp) function expm1_dp(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: term, sumv
      integer :: k
      if (abs(x) > 1.0e-5_dp) then
         v = exp(x) - 1.0_dp
      else
         term = x
         sumv = x
         do k = 2, 30
            term = term*x/real(k,dp)
            sumv = sumv + term
            if (abs(term) <= epsilon(1.0_dp)*max(1.0_dp,abs(sumv))) exit
         end do
         v = sumv
      end if
   end function expm1_dp

   pure elemental function d_0(log_p) result(v)
      logical, intent(in) :: log_p
      real(dp) :: v
      if (log_p) then
         v = ieee_value(0.0_dp, ieee_negative_inf)
      else
         v = 0.0_dp
      end if
   end function d_0

   pure elemental function d_1(log_p) result(v)
      logical, intent(in) :: log_p
      real(dp) :: v
      if (log_p) then
         v = 0.0_dp
      else
         v = 1.0_dp
      end if
   end function d_1

   pure elemental function dt_0(lower_tail, log_p) result(v)
      logical, intent(in) :: lower_tail, log_p
      real(dp) :: v
      if (lower_tail) then
         v = d_0(log_p)
      else
         v = d_1(log_p)
      end if
   end function dt_0

   pure elemental function dt_1(lower_tail, log_p) result(v)
      logical, intent(in) :: lower_tail, log_p
      real(dp) :: v
      if (lower_tail) then
         v = d_1(log_p)
      else
         v = d_0(log_p)
      end if
   end function dt_1

   pure elemental real(dp) function d_lval(p, lower_tail) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail
      if (lower_tail) then
         v = p
      else
         v = 1.0_dp - p
      end if
   end function d_lval

   pure elemental real(dp) function d_cval(p, lower_tail) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail
      if (lower_tail) then
         v = 1.0_dp - p
      else
         v = p
      end if
   end function d_cval

   pure elemental real(dp) function d_val(x, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in) :: log_p
      if (log_p) then
         v = log(x)
      else
         v = x
      end if
   end function d_val

   pure elemental real(dp) function d_qiv(p, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: log_p
      if (log_p) then
         v = exp(p)
      else
         v = p
      end if
   end function d_qiv

   pure elemental real(dp) function d_exp(x, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in) :: log_p
      if (log_p) then
         v = x
      else
         v = exp(x)
      end if
   end function d_exp

   pure elemental real(dp) function d_log(p, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: log_p
      if (log_p) then
         v = p
      else
         v = log(p)
      end if
   end function d_log

   pure elemental real(dp) function d_clog(p, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: log_p
      if (log_p) then
         v = log1p_dp(-exp(p))
      else
         v = (0.5_dp - p) + 0.5_dp
      end if
   end function d_clog

   pure elemental real(dp) function log1mexp(x) result(v)
      real(dp), intent(in) :: x
      if (x < 0.0_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x <= m_ln2) then
         v = log(-expm1_dp(-x))
      else
         v = log1p_dp(-exp(-x))
      end if
   end function log1mexp

   pure elemental real(dp) function log1pexp(x) result(v)
      real(dp), intent(in) :: x
      if (x <= -37.0_dp) then
         v = exp(x)
      else if (x <= 18.0_dp) then
         v = log1p_dp(exp(x))
      else if (x <= 33.3_dp) then
         v = x + exp(-x)
      else
         v = x
      end if
   end function log1pexp

   pure elemental real(dp) function d_lexp(x, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in) :: log_p
      if (log_p) then
         v = log1mexp(-x)
      else
         v = log1p_dp(-x)
      end if
   end function d_lexp

   pure elemental real(dp) function dt_val(x, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in) :: lower_tail, log_p
      v = d_val(d_lval(x, lower_tail), log_p)
   end function dt_val

   pure elemental real(dp) function dt_cval(x, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in) :: lower_tail, log_p
      v = d_val(d_cval(x, lower_tail), log_p)
   end function dt_cval

   pure elemental real(dp) function dt_qiv(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail, log_p
      if (log_p) then
         if (lower_tail) then
            v = exp(p)
         else
            v = -expm1_dp(p)
         end if
      else
         v = d_lval(p, lower_tail)
      end if
   end function dt_qiv

   pure elemental real(dp) function dt_civ(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail, log_p
      if (log_p) then
         if (lower_tail) then
            v = -expm1_dp(p)
         else
            v = exp(p)
         end if
      else
         v = d_cval(p, lower_tail)
      end if
   end function dt_civ

   pure elemental real(dp) function dt_log(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail, log_p
      if (lower_tail) then
         v = d_log(p, log_p)
      else
         v = d_lexp(p, log_p)
      end if
   end function dt_log

   pure elemental real(dp) function dt_clog(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail, log_p
      if (lower_tail) then
         v = d_lexp(p, log_p)
      else
         v = d_log(p, log_p)
      end if
   end function dt_clog

   pure elemental real(dp) function dt_log_known(p, lower_tail) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail
      if (lower_tail) then
         v = p
      else
         v = log1mexp(-p)
      end if
   end function dt_log_known

   pure elemental real(dp) function clamp01(x) result(v)
      real(dp), intent(in) :: x
      v = max(0.0_dp, min(1.0_dp, x))
   end function clamp01

   pure elemental real(dp) function prob_from_input(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail, log_p
      v = dt_qiv(p, lower_tail, log_p)
   end function prob_from_input

   pure elemental real(dp) function lower_prob_log(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in) :: lower_tail, log_p
      v = dt_log(p, lower_tail, log_p)
   end function lower_prob_log

   pure elemental real(dp) function prob_output(p_lower, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p_lower
      logical, intent(in) :: lower_tail, log_p
      real(dp) :: p
      p = clamp01(p_lower)
      if (lower_tail) then
         if (log_p) then
            if (p == 0.0_dp) then
               v = ieee_value(0.0_dp, ieee_negative_inf)
            else
               v = log(p)
            end if
         else
            v = p
         end if
      else
         if (log_p) then
            if (p == 1.0_dp) then
               v = ieee_value(0.0_dp, ieee_negative_inf)
            else
               v = log1p_dp(-p)
            end if
         else
            v = 1.0_dp - p
         end if
      end if
   end function prob_output

   pure elemental real(dp) function logspace_add(logx, logy) result(v)
      real(dp), intent(in) :: logx, logy
      real(dp) :: m
      if (logx == ieee_value(0.0_dp, ieee_negative_inf)) then
         v = logy
      else if (logy == ieee_value(0.0_dp, ieee_negative_inf)) then
         v = logx
      else
         m = max(logx, logy)
         v = m + log1p_dp(exp(min(logx, logy) - m))
      end if
   end function logspace_add

   pure elemental real(dp) function logspace_sub(logx, logy) result(v)
      real(dp), intent(in) :: logx, logy
      if (logy > logx) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (logy == ieee_value(0.0_dp, ieee_negative_inf)) then
         v = logx
      else if (logy == logx) then
         v = ieee_value(0.0_dp, ieee_negative_inf)
      else
         v = logx + log1p_dp(-exp(logy-logx))
      end if
   end function logspace_sub

   pure elemental real(dp) function pow1p(x, y) result(v)
      real(dp), intent(in) :: x, y
      integer :: iy
      if (y >= 0.0_dp .and. y <= 4.0_dp .and. y == anint(y)) then
         iy = int(y)
         select case (iy)
         case (0)
            v = 1.0_dp
         case (1)
            v = 1.0_dp + x
         case (2)
            v = x*(x + 2.0_dp) + 1.0_dp
         case (3)
            v = x*(x*(x + 3.0_dp) + 3.0_dp) + 1.0_dp
         case (4)
            v = x*(x*(x*(x + 4.0_dp) + 6.0_dp) + 4.0_dp) + 1.0_dp
         end select
      else if (((1.0_dp + x) - 1.0_dp == x) .or. abs(x) > 0.5_dp) then
         v = (1.0_dp + x)**y
      else
         v = exp(y*log1p_dp(x))
      end if
   end function pow1p

   pure elemental real(dp) function dpq_pow(x, y) result(v)
      real(dp), intent(in) :: x, y
      if (x == 1.0_dp .or. y == 0.0_dp) then
         v = 1.0_dp
      else
         v = x**y
      end if
   end function dpq_pow

   pure elemental real(dp) function pow_di(x, n) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: base
      integer :: k
      v = 1.0_dp
      if (n == 0) return
      base = x
      k = abs(n)
      do while (k > 0)
         if (iand(k,1) /= 0) v = v*base
         k = ishft(k,-1)
         if (k > 0) base = base*base
      end do
      if (n < 0) v = 1.0_dp/v
   end function pow_di

   pure elemental real(dp) function expm1x_tser(x, k) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: k
      real(dp) :: p
      integer :: j
      if (k < 1) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      p = x / real(k+1,dp)
      do j = k, 2, -1
         p = (1.0_dp + p)*x/real(j,dp)
      end do
      v = x*p
   end function expm1x_tser

   pure elemental real(dp) function expm1x(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: ax
      ax = abs(x)
      if (ax < 4.4e-8_dp) then
         v = expm1x_tser(x,2)
      else if (ax < 0.10_dp) then
         v = expm1x_tser(x,9)
      else if (ax < 0.385_dp) then
         v = expm1x_tser(x,12)
      else if (ax < 1.1_dp) then
         v = expm1x_tser(x,17)
      else if (ax < 2.0_dp) then
         v = expm1_dp(x) - x
      else
         v = exp(x) - 1.0_dp - x
      end if
   end function expm1x


   pure elemental real(dp) function rexpm1(x) result(v)
      real(dp), intent(in) :: x
      real(dp), parameter :: p1=9.14041914819518e-10_dp, p2=0.0238082361044469_dp
      real(dp), parameter :: q1=-0.499999999085958_dp, q2=0.107141568980644_dp
      real(dp), parameter :: q3=-0.0119041179760821_dp, q4=5.95130811860248e-4_dp
      real(dp) :: w
      if (abs(x) <= 0.15_dp) then
         v=x*(((p2*x+p1)*x+1.0_dp)/((((q4*x+q3)*x+q2)*x+q1)*x+1.0_dp))
      else
         w=exp(x)
         if(x>0.0_dp)then
            v=w*(0.5_dp-1.0_dp/w+0.5_dp)
         else
            v=(w-0.5_dp)-0.5_dp
         end if
      end if
   end function rexpm1

   pure elemental real(dp) function rlog1(x) result(v)
      real(dp), intent(in) :: x
      real(dp), parameter :: p0=0.333333333333333_dp, p1=-0.224696413112536_dp
      real(dp), parameter :: p2=0.00620886815375787_dp, q1=-1.27408923933623_dp
      real(dp), parameter :: q2=0.354508718369557_dp
      real(dp) :: h,w1,r,t,w
      if(x>=-0.39_dp .and. x<0.57_dp)then
         h=x
         w1=0.0_dp
         if(x< -0.18_dp)then
            h=(x+0.3_dp)/0.7_dp
            w1=0.0566749439387324_dp-0.3_dp*h
         else if(x>0.18_dp)then
            h=0.75_dp*x-0.25_dp
            w1=0.0456512608815524_dp+h/3.0_dp
         end if
         r=h/(h+2.0_dp)
         t=r*r
         w=((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0_dp)
         v=2.0_dp*t*(1.0_dp/(1.0_dp-r)-r*w)+w1
      else
         v=x-log((x+0.5_dp)+0.5_dp)
      end if
   end function rlog1

   pure elemental real(dp) function log4p1p(p,log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: log_p
      logical :: lp
      lp=.false.
      if(present(log_p))lp=log_p
      if(lp)then
         v=log(-4.0_dp*expm1_dp(p))+p
      else
         v=log(4.0_dp*p*(1.0_dp-p))
      end if
   end function log4p1p

   pure elemental real(dp) function p1l1p(t) result(v)
      real(dp), intent(in) :: t
      v = log1pmx(t) + t*log1p_dp(t)
   end function p1l1p

   pure elemental real(dp) function p1l1ser(t, k) result(v)
      real(dp), intent(in) :: t
      integer, intent(in) :: k
      real(dp) :: term
      integer :: j
      if (t <= -1.0_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      ! (1+t)log(1+t)-t = sum_{j=2} inf (-1)^j t^j/[j(j-1)]
      v = 0.0_dp
      term = t*t
      do j = 2, max(2,k+1)
         if (mod(j,2) == 0) then
            v = v + term/(real(j,dp)*real(j-1,dp))
         else
            v = v - term/(real(j,dp)*real(j-1,dp))
         end if
         term = term*t
      end do
   end function p1l1ser

   pure elemental real(dp) function p1l1(t) result(v)
      real(dp), intent(in) :: t
      if (t <= -1.0_dp) then
         if (t == -1.0_dp) then
            v = 1.0_dp
         else
            v = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      else if (abs(t) < 0.25_dp) then
         v = p1l1p(t)
      else
         v = (1.0_dp+t)*log1p_dp(t)-t
      end if
   end function p1l1

   pure elemental real(dp) function logr(x, a) result(v)
      real(dp), intent(in) :: x, a
      if (x < 0.0_dp .or. x+a <= 0.0_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (a == 0.0_dp) then
         v = 0.0_dp
      else if (x == 0.0_dp) then
         v = ieee_value(0.0_dp, ieee_negative_inf)
      else
         v = -log1p_dp(a/x)
      end if
   end function logr

   pure real(dp) function logcf(x, i, d, eps, maxit) result(v)
      real(dp), intent(in) :: x, eps
      integer, intent(in) :: i, d
      integer, intent(in), optional :: maxit
      real(dp) :: c1, c2, c4, a1, b1, a2, b2, c3
      integer :: k, imax
      imax = 10000
      if (present(maxit)) imax = maxit
      if (x == 0.0_dp) then
         v = 1.0_dp/real(i,dp)
         return
      end if
      c1 = 2.0_dp*real(d,dp)
      c2 = real(i,dp) + real(d,dp)
      c4 = c2 + real(d,dp)
      a1 = c2
      b1 = real(i,dp)*(c2-real(i,dp)*x)
      b2 = real(d,dp)*real(d,dp)*x
      a2 = c4*c2-b2
      b2 = c4*b1-real(i,dp)*b2
      if (b2 /= 0.0_dp) then
         c3 = a2/b2
      else
         c3 = 0.0_dp
      end if
      do k = 1, imax
         c3 = c3*c1
         c1 = c1 + real(d,dp)
         c2 = c2 + real(d,dp)
         c4 = c4 + real(d,dp)
         a1 = c4*a2-c3*a1
         b1 = c4*b2-c3*b1
         c3 = c1*c2*x
         a2 = c4*a1-c3*a2
         b2 = c4*b1-c3*b2
         if (b2 /= 0.0_dp) then
            if (abs(a2/b2-a1/b1) <= abs(eps*a2/b2)) exit
         end if
      end do
      v = a2/b2
   end function logcf

   pure elemental real(dp) function log1pmx(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: term, y
      if (x < -1.0_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == -1.0_dp) then
         v = ieee_value(0.0_dp, ieee_negative_inf)
      else if (x >= 2.0_dp**59) then
         v = -x
      else if (x > 1.0_dp .or. x < -0.79149064_dp) then
         v = log1p_dp(x)-x
      else
         term = x/(2.0_dp+x)
         y = term*term
         if (abs(x) <= 0.01_dp) then
            v = term*((((2.0_dp/9.0_dp*y + 2.0_dp/7.0_dp)*y + &
               2.0_dp/5.0_dp)*y + 2.0_dp/3.0_dp)*y - x)
         else
            v = term*(2.0_dp*y*logcf(y,3,2,1.0e-14_dp)-x)
         end if
      end if
   end function log1pmx

   pure elemental real(dp) function lgamma1p(a) result(v)
      real(dp), intent(in) :: a
      real(dp), parameter :: euler = &
         0.577215664901532860606512090082402431_dp
      real(dp), parameter :: coeffs(40) = [ &
         0.3224670334241132182362075833230126_dp, &
         0.0673523010531980951332460538371500_dp, &
         0.02058080842778454787900092413529198_dp, &
         0.007385551028673985266273097291406834_dp, &
         0.002890510330741523285752988298486755_dp, &
         0.001192753911703260977113935692828109_dp, &
         0.0005096695247430424223356548135815582_dp, &
         0.0002231547584535793797614188036013401_dp, &
         0.00009945751278180853371459589003190170_dp, &
         0.00004492623673813314170020750240635786_dp, &
         0.00002050721277567069155316650397830591_dp, &
         0.000009439488275268395903987425104415055_dp, &
         0.000004374866789907487804181793223952411_dp, &
         0.000002039215753801366236781900709670839_dp, &
         0.0000009551412130407419832857179772951265_dp, &
         0.0000004492469198764566043294290331193655_dp, &
         0.0000002120718480555466586923135901077628_dp, &
         0.0000001004322482396809960872083050053344_dp, &
         4.769810169363980565760193417246730e-8_dp, &
         2.271109460894316491031998116062124e-8_dp, &
         1.083865921489695409107491757968159e-8_dp, &
         5.183475041970046655121248647057669e-9_dp, &
         2.483674543802478317185008663991718e-9_dp, &
         1.192140140586091207442548202774640e-9_dp, &
         5.731367241678862013330194857961011e-10_dp, &
         2.759522885124233145178149692816341e-10_dp, &
         1.330476437424448948149715720858008e-10_dp, &
         6.422964563838100022082448087640e-11_dp, &
         3.104424774732227276239215783401e-11_dp, &
         1.502138408075414217093301048781e-11_dp, &
         7.275974480239079662504549924814e-12_dp, &
         3.527742476575915083615075225654e-12_dp, &
         1.711991790559617908601084114443e-12_dp, &
         8.315385841420284819798357793954e-13_dp, &
         4.042200525289440065536008957033e-13_dp, &
         1.966475631096616490411045679010e-13_dp, &
         9.573630387838555763782200936509e-14_dp, &
         4.664076026428374224576492564978e-14_dp, &
         2.273736960065972320633279596737e-14_dp, &
         1.109139947083452201658320007192e-14_dp ]
      real(dp), parameter :: c = 2.273736845824652515226821577978691e-13_dp
      real(dp) :: lgam
      integer :: i
      if (abs(a) >= 0.5_dp) then
         v = r_lgamma(1.0_dp+a)
      else
         lgam = c*logcf(-a/2.0_dp,42,1,1.0e-14_dp)
         do i = 40, 1, -1
            lgam = coeffs(i)-a*lgam
         end do
         v = (a*lgam-euler)*a-log1pmx(a)
      end if
   end function lgamma1p

   pure elemental real(dp) function lgamma1p_series(x, k) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: k
      real(dp), parameter :: euler = &
         0.577215664901532860606512090082402431_dp
      real(dp), parameter :: zeta3 = 1.20205690315959428539973816151145_dp
      real(dp), parameter :: zeta5 = 1.03692775514336992633136548645703_dp
      real(dp), parameter :: zeta7 = 1.00834927738192282683979754984980_dp
      real(dp), parameter :: zeta9 = 1.00200839282608221441785276923241_dp
      real(dp) :: xx
      integer :: kk
      ! First 10 coefficients of log Gamma(1+x): (-1)^n zeta(n)/n, n>=2.
      v = -euler*x
      xx = x*x
      kk = min(max(k,1),10)
      if (kk >= 2) v = v + (acos(-1.0_dp)**2/12.0_dp)*xx
      if (kk >= 3) v = v - (zeta3/3.0_dp)*xx*x
      if (kk >= 4) v = v + (acos(-1.0_dp)**4/360.0_dp)*xx*xx
      if (kk >= 5) v = v - (zeta5/5.0_dp)*xx*xx*x
      if (kk >= 6) v = v + (acos(-1.0_dp)**6/5670.0_dp)*xx**3
      if (kk >= 7) v = v - (zeta7/7.0_dp)*xx**3*x
      if (kk >= 8) v = v + (acos(-1.0_dp)**8/75600.0_dp)*xx**4
      if (kk >= 9) v = v - (zeta9/9.0_dp)*xx**4*x
      if (kk >= 10) v = v + (acos(-1.0_dp)**10/935550.0_dp)*xx**5
   end function lgamma1p_series

   pure function lsum(lx) result(v)
      real(dp), intent(in) :: lx(:)
      real(dp) :: v, m
      if (size(lx) == 0) then
         v = ieee_value(0.0_dp, ieee_negative_inf)
         return
      end if
      m = maxval(lx)
      if (m == ieee_value(0.0_dp, ieee_positive_inf)) then
         v = m
      else if (.not. ieee_is_finite(m)) then
         v = m
      else
         v = m + log(sum(exp(lx-m)))
      end if
   end function lsum

   pure function lssum(lxabs, signs) result(v)
      real(dp), intent(in) :: lxabs(:)
      integer, intent(in) :: signs(:)
      real(dp) :: v, mp, mn, sp, sn
      if (size(lxabs) /= size(signs)) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (count(signs > 0) == 0) then
         if (count(signs < 0) == 0) then
            v = ieee_value(0.0_dp, ieee_negative_inf)
         else
            mn = maxval(lxabs, mask=signs<0)
            sn = sum(exp(lxabs-mn), mask=signs<0)
            v = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
         return
      end if
      mp = maxval(lxabs, mask=signs>0)
      sp = sum(exp(lxabs-mp), mask=signs>0)
      if (count(signs < 0) == 0) then
         v = mp + log(sp)
         return
      end if
      mn = maxval(lxabs, mask=signs<0)
      sn = sum(exp(lxabs-mn), mask=signs<0)
      if (mp >= mn) then
         sn = exp(mn-mp)*sn
         if (sp <= sn) then
            v = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            v = mp + log(sp-sn)
         end if
      else
         sp = exp(mp-mn)*sp
         if (sp <= sn) then
            v = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            v = mn + log(sp-sn)
         end if
      end if
   end function lssum

   pure function chebyshev_eval(x, a, n) result(v)
      real(dp), intent(in) :: x, a(:)
      integer, intent(in), optional :: n
      real(dp) :: v, b0, b1, b2, twox
      integer :: i, nn
      nn = size(a)
      if (present(n)) nn = min(nn,n)
      if (nn <= 0 .or. abs(x) > 1.1_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      twox = 2.0_dp*x
      b0 = 0.0_dp
      b1 = 0.0_dp
      b2 = 0.0_dp
      do i = nn, 1, -1
         b2 = b1
         b1 = b0
         b0 = twox*b1-b2+a(i)
      end do
      v = 0.5_dp*(b0-b2)
   end function chebyshev_eval

   pure function chebyshev_poly(x, n) result(t)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: t(0:n)
      integer :: j
      t(0) = 1.0_dp
      if (n >= 1) t(1) = x
      do j = 2, n
         t(j) = 2.0_dp*x*t(j-1)-t(j-2)
      end do
   end function chebyshev_poly

   pure function chebyshev_nc(x, n) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: v
      real(dp) :: t(0:max(0,n))
      t = chebyshev_poly(x,max(0,n))
      v = t(max(0,n))
   end function chebyshev_nc

   pure subroutine dpq_frexp(x, frac, expo)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: frac
      integer, intent(out) :: expo
      if (x == 0.0_dp) then
         frac = 0.0_dp
         expo = 0
      else
         frac = fraction(x)
         expo = exponent(x)
      end if
   end subroutine dpq_frexp

   pure elemental real(dp) function dpq_ldexp(x, expo) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: expo
      v = scale(x,expo)
   end function dpq_ldexp

   pure subroutine dpq_modf(x, ipart, fpart)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: ipart, fpart
      ipart = aint(x)
      fpart = x-ipart
   end subroutine dpq_modf

   pure elemental real(dp) function dpsifn_scalar(x, deriv) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: deriv
      v = r_psigamma(x,deriv)
   end function dpsifn_scalar

end module dpq_core
