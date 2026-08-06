! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_math
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use sn_kinds, only : dp, pi, sqrt_two, sqrt_two_pi, log_two_pi, tiny_dp, huge_dp
  implicit none
  private

  public :: normal_pdf, normal_logpdf, normal_cdf, normal_logcdf, normal_quantile
  public :: student_t_pdf, student_t_logpdf, student_t_cdf, student_t_quantile
  public :: chi_square_cdf, chi_square_quantile, f_cdf, f_quantile
  public :: regularized_beta, regularized_gamma_p
  public :: adaptive_simpson, clamp_probability, log1mexp
  public :: log_sum_exp, finite_real, log1p_dp, expm1_dp

  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
  end interface

contains

  pure real(dp) function log1p_dp(x) result(value)
    real(dp), intent(in) :: x
    if (abs(x) > 1.0e-4_dp) then
      value = log(1.0_dp+x)
    else
      value = x*(1.0_dp+x*(-0.5_dp+x*(1.0_dp/3.0_dp+x*(-0.25_dp+x*0.2_dp))))
    end if
  end function log1p_dp

  pure real(dp) function expm1_dp(x) result(value)
    real(dp), intent(in) :: x
    if (abs(x) > 1.0e-5_dp) then
      value = exp(x)-1.0_dp
    else
      value = x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x*(1.0_dp/24.0_dp+x/120.0_dp))))
    end if
  end function expm1_dp

  pure real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp*x*x)/sqrt_two_pi
  end function normal_pdf

  pure real(dp) function normal_logpdf(x) result(value)
    real(dp), intent(in) :: x
    value = -0.5_dp*(log_two_pi + x*x)
  end function normal_logpdf

  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp*erfc(-x/sqrt_two)
  end function normal_cdf

  pure real(dp) function normal_logcdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: z, corr

    if (x > -8.0_dp) then
      value = log(max(normal_cdf(x), tiny_dp))
    else
      z = -x
      corr = 1.0_dp - 1.0_dp/(z*z) + 3.0_dp/(z**4) - 15.0_dp/(z**6)
      corr = max(corr, tiny_dp)
      value = -0.5_dp*z*z - log(z) - 0.5_dp*log_two_pi + log(corr)
    end if
  end function normal_logcdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
      -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
      -3.066479806614716e+01_dp, 2.506628277459239e+00_dp]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
      -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
      -1.328068155288572e+01_dp]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
      -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
       4.374664141464968e+00_dp, 2.938163982698783e+00_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp]
    real(dp), parameter :: p_low = 0.02425_dp
    real(dp), parameter :: p_high = 1.0_dp - p_low
    real(dp) :: q, r, e, u

    if (p <= 0.0_dp) then
      x = -huge_dp
      return
    else if (p >= 1.0_dp) then
      x = huge_dp
      return
    end if

    if (p < p_low) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= p_high) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if

    ! One Halley correction provides near machine precision.
    e = normal_cdf(x) - p
    u = e/normal_pdf(x)
    x = x - u/(1.0_dp + 0.5_dp*x*u)
  end function normal_quantile

  pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: max_iter = 400
    real(dp), parameter :: eps = 4.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny_dp/eps
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
    do m = 1, max_iter
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x/ &
           ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
           ((a+real(m2,dp))*(qap+real(m2,dp)))
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
  end function beta_continued_fraction

  pure real(dp) function regularized_beta(a, b, x) result(value)
    real(dp), intent(in) :: a, b, x
    real(dp) :: bt

    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      value = 1.0_dp
      return
    end if

    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log1p_dp(-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      value = bt*beta_continued_fraction(a,b,x)/a
    else
      value = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
    end if
    value = min(1.0_dp, max(0.0_dp, value))
  end function regularized_beta

  pure real(dp) function regularized_gamma_p(a, x) result(value)
    real(dp), intent(in) :: a, x
    integer, parameter :: max_iter = 500
    real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny_dp/eps
    integer :: n
    real(dp) :: ap, del, sum, b, c, d, h, an

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (abs(x) <= tiny_dp) then
      value = 0.0_dp
      return
    end if

    if (x < a+1.0_dp) then
      ap = a
      sum = 1.0_dp/a
      del = sum
      do n = 1, max_iter
        ap = ap + 1.0_dp
        del = del*x/ap
        sum = sum + del
        if (abs(del) <= abs(sum)*eps) exit
      end do
      value = sum*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n = 1, max_iter
        an = -real(n,dp)*(real(n,dp)-a)
        b = b + 2.0_dp
        d = an*d+b
        if (abs(d) < fpmin) d = fpmin
        c = b+an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      value = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    value = min(1.0_dp, max(0.0_dp, value))
  end function regularized_gamma_p

  pure real(dp) function student_t_logpdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    if (nu <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
    else
      value = log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu) &
              -0.5_dp*log(nu*pi)-0.5_dp*(nu+1.0_dp)*log1p_dp(x*x/nu)
    end if
  end function student_t_logpdf

  pure real(dp) function student_t_pdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    value = exp(student_t_logpdf(x,nu))
  end function student_t_pdf

  pure real(dp) function student_t_cdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    real(dp) :: ib, z
    if (nu <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (abs(x) <= tiny_dp) then
      value = 0.5_dp
      return
    end if
    z = nu/(nu+x*x)
    ib = regularized_beta(0.5_dp*nu,0.5_dp,z)
    if (x > 0.0_dp) then
      value = 1.0_dp-0.5_dp*ib
    else
      value = 0.5_dp*ib
    end if
  end function student_t_cdf

  real(dp) function student_t_quantile(p, nu) result(x)
    real(dp), intent(in) :: p, nu
    real(dp) :: lo, hi, mid, pmid
    integer :: iter

    if (p <= 0.0_dp) then
      x = -huge_dp
      return
    else if (p >= 1.0_dp) then
      x = huge_dp
      return
    else if (nu <= 0.0_dp) then
      x = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (abs(p-0.5_dp) <= epsilon(1.0_dp)) then
      x = 0.0_dp
      return
    end if

    lo = -1.0_dp
    hi = 1.0_dp
    do while (student_t_cdf(lo,nu) > p)
      lo = 2.0_dp*lo
      if (lo < -1.0e12_dp) exit
    end do
    do while (student_t_cdf(hi,nu) < p)
      hi = 2.0_dp*hi
      if (hi > 1.0e12_dp) exit
    end do
    do iter = 1, 160
      mid = 0.5_dp*(lo+hi)
      pmid = student_t_cdf(mid,nu)
      if (pmid < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) <= 2.0e-13_dp*max(1.0_dp,abs(mid))) exit
    end do
    x = 0.5_dp*(lo+hi)
  end function student_t_quantile

  pure real(dp) function chi_square_cdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    if (x <= 0.0_dp) then
      value = 0.0_dp
    else
      value = regularized_gamma_p(0.5_dp*nu,0.5_dp*x)
    end if
  end function chi_square_cdf

  real(dp) function chi_square_quantile(p, nu) result(x)
    real(dp), intent(in) :: p, nu
    real(dp) :: lo, hi, mid
    integer :: iter
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge_dp
      return
    end if
    lo = 0.0_dp
    hi = max(1.0_dp,nu)
    do while (chi_square_cdf(hi,nu) < p)
      hi = 2.0_dp*hi
    end do
    do iter=1,160
      mid = 0.5_dp*(lo+hi)
      if (chi_square_cdf(mid,nu) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function chi_square_quantile

  pure real(dp) function f_cdf(x, d1, d2) result(value)
    real(dp), intent(in) :: x, d1, d2
    if (x <= 0.0_dp) then
      value = 0.0_dp
    else
      value = regularized_beta(0.5_dp*d1,0.5_dp*d2,d1*x/(d1*x+d2))
    end if
  end function f_cdf

  real(dp) function f_quantile(p, d1, d2) result(x)
    real(dp), intent(in) :: p, d1, d2
    real(dp) :: lo, hi, mid
    integer :: iter
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge_dp
      return
    end if
    lo = 0.0_dp
    hi = 1.0_dp
    do while (f_cdf(hi,d1,d2) < p)
      hi = 2.0_dp*hi
    end do
    do iter=1,160
      mid = 0.5_dp*(lo+hi)
      if (f_cdf(mid,d1,d2) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function f_quantile

  real(dp) function adaptive_simpson(f, a, b, tol, max_depth) result(value)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_depth
    real(dp) :: eps, fa, fb, fm, whole
    integer :: depth

    eps = 1.0e-10_dp
    if (present(tol)) eps = max(tol,10.0_dp*epsilon(1.0_dp))
    depth = 20
    if (present(max_depth)) depth = max_depth
    fa = f(a)
    fb = f(b)
    fm = f(0.5_dp*(a+b))
    whole = (b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    value = recurse(a,b,fa,fb,fm,whole,eps,depth)
  contains
    recursive real(dp) function recurse(x0,x1,f0,f1,fm0,s,eps0,dleft) result(ans)
      real(dp), intent(in) :: x0,x1,f0,f1,fm0,s,eps0
      integer, intent(in) :: dleft
      real(dp) :: xm, xl, xr, fl, fr, sl, sr, delta
      xm = 0.5_dp*(x0+x1)
      xl = 0.5_dp*(x0+xm)
      xr = 0.5_dp*(xm+x1)
      fl = f(xl)
      fr = f(xr)
      sl = (xm-x0)*(f0+4.0_dp*fl+fm0)/6.0_dp
      sr = (x1-xm)*(fm0+4.0_dp*fr+f1)/6.0_dp
      delta = sl+sr-s
      if (dleft <= 0 .or. abs(delta) <= 15.0_dp*eps0) then
        ans = sl+sr+delta/15.0_dp
      else
        ans = recurse(x0,xm,f0,fm0,fl,sl,0.5_dp*eps0,dleft-1) + &
              recurse(xm,x1,fm0,f1,fr,sr,0.5_dp*eps0,dleft-1)
      end if
    end function recurse
  end function adaptive_simpson

  pure real(dp) function clamp_probability(p) result(value)
    real(dp), intent(in) :: p
    value = min(1.0_dp,max(0.0_dp,p))
  end function clamp_probability

  pure real(dp) function log1mexp(x) result(value)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
    else if (x < -log(2.0_dp)) then
      value = log1p_dp(-exp(x))
    else
      value = log(-expm1_dp(x))
    end if
  end function log1mexp

  pure real(dp) function log_sum_exp(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: xmax
    if (size(x) == 0) then
      value = -huge_dp
      return
    end if
    xmax = maxval(x)
    if (.not. finite_real(xmax)) then
      value = xmax
    else
      value = xmax + log(sum(exp(x-xmax)))
    end if
  end function log_sum_exp

  pure logical function finite_real(x) result(ok)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp), intent(in) :: x
    ok = ieee_is_finite(x)
  end function finite_real

end module sn_math
