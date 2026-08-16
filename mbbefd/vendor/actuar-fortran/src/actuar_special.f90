! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_special
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use actuar_kinds, only : dp, pi, sqrt2, huge_dp
  implicit none
  private
  public :: nan_dp, clamp01, log1pexp, log1mexp
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: regularized_gamma_p, regularized_gamma_q, gamma_quantile
  public :: regularized_beta, beta_quantile, log_beta_fn
  public :: bessel_k_integral, adaptive_simpson

  abstract interface
    pure function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface

contains

  pure function nan_dp() result(x)
    real(dp) :: x
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure function clamp01(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = max(0.0_dp, min(1.0_dp, x))
  end function clamp01

  pure function log1pexp(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x > 35.0_dp) then
      y = x + exp(-x)
    else if (x < -35.0_dp) then
      y = exp(x)
    else
      y = log(1.0_dp + exp(x))
    end if
  end function log1pexp

  pure function log1mexp(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x >= 0.0_dp) then
      y = nan_dp()
    else if (x < -log(2.0_dp)) then
      y = log(1.0_dp - exp(x))
    else
      y = log(1.0_dp-exp(x))
    end if
  end function log1mexp

  pure function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  pure function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 0.5_dp * erfc(-x/sqrt2)
  end function normal_cdf

  pure function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp
    real(dp), parameter :: phigh = 1.0_dp - plow

    if (p <= 0.0_dp) then
      x = -huge_dp
      return
    else if (p >= 1.0_dp) then
      x = huge_dp
      return
    end if

    if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if

    x = x - (normal_cdf(x)-p) / max(normal_pdf(x), tiny(1.0_dp))
  end function normal_quantile

  pure function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p, sumv, del, ap, gln, b, c, d, h, an
    integer :: n
    real(dp), parameter :: eps = 4.0e-15_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = nan_dp()
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, 10000
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/max(b, fpmin)
      h = d
      do n = 1, 10000
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
      p = 1.0_dp - exp(-x+a*log(x)-gln)*h
    end if
    p = clamp01(p)
  end function regularized_gamma_p

  pure function regularized_gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    real(dp) :: q
    q = 1.0_dp - regularized_gamma_p(a, x)
    q = clamp01(q)
  end function regularized_gamma_q

  function gamma_quantile(p, shape) result(x)
    real(dp), intent(in) :: p, shape
    real(dp) :: x, lo, hi, mid
    integer :: iter
    if (shape <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp(); return
    end if
    if (p <= 0.0_dp) then
      x = 0.0_dp; return
    else if (p >= 1.0_dp) then
      x = huge_dp; return
    end if
    lo = 0.0_dp
    hi = max(1.0_dp, shape + 8.0_dp*sqrt(shape) + 20.0_dp)
    do while (regularized_gamma_p(shape, hi) < p)
      hi = 2.0_dp*hi
      if (hi > huge_dp/4.0_dp) exit
    end do
    do iter = 1, 160
      mid = 0.5_dp*(lo+hi)
      if (regularized_gamma_p(shape, mid) < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) <= 4.0e-14_dp*max(1.0_dp,mid)) exit
    end do
    x = 0.5_dp*(lo+hi)
  end function gamma_quantile

  pure function log_beta_fn(a, b) result(x)
    real(dp), intent(in) :: a, b
    real(dp) :: x
    x = log_gamma(a) + log_gamma(b) - log_gamma(a+b)
  end function log_beta_fn

  pure function beta_cf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    real(dp) :: cf, qab, qap, qam, c, d, h, aa, del
    integer :: m, m2
    real(dp), parameter :: eps = 4.0e-15_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    qab = a+b
    qap = a+1.0_dp
    qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, 10000
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
           ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
           ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_cf

  pure function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: p, bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = nan_dp(); return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp; return
    else if (x >= 1.0_dp) then
      p = 1.0_dp; return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      p = bt*beta_cf(a,b,x)/a
    else
      p = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    p = clamp01(p)
  end function regularized_beta

  function beta_quantile(p, a, b) result(x)
    real(dp), intent(in) :: p, a, b
    real(dp) :: x, lo, hi, mid
    integer :: iter
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp(); return
    end if
    if (p <= 0.0_dp) then
      x = 0.0_dp; return
    else if (p >= 1.0_dp) then
      x = 1.0_dp; return
    end if
    lo = 0.0_dp; hi = 1.0_dp
    do iter = 1, 180
      mid = 0.5_dp*(lo+hi)
      if (regularized_beta(mid,a,b) < p) then
        lo = mid
      else
        hi = mid
      end if
      if (hi-lo < 2.0e-14_dp) exit
    end do
    x = 0.5_dp*(lo+hi)
  end function beta_quantile

  function adaptive_simpson(f, a, b, tol, max_depth) result(value)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a, b, tol
    integer, intent(in), optional :: max_depth
    real(dp) :: value, fa, fb, fm, whole
    integer :: depth
    depth = 20
    if (present(max_depth)) depth = max_depth
    fa = f(a); fb = f(b); fm = f(0.5_dp*(a+b))
    whole = (b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    value = simpson_rec(f,a,b,fa,fm,fb,whole,tol,depth)
  end function adaptive_simpson

  recursive function simpson_rec(f,a,b,fa,fm,fb,whole,tol,depth) result(v)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b,fa,fm,fb,whole,tol
    integer, intent(in) :: depth
    real(dp) :: v,c,lm,rm,flm,frm,left,right,delta
    c = 0.5_dp*(a+b)
    lm = 0.5_dp*(a+c); rm = 0.5_dp*(c+b)
    flm = f(lm); frm = f(rm)
    left = (c-a)*(fa+4.0_dp*flm+fm)/6.0_dp
    right = (b-c)*(fm+4.0_dp*frm+fb)/6.0_dp
    delta = left+right-whole
    if (depth <= 0 .or. abs(delta) <= 15.0_dp*tol) then
      v = left+right+delta/15.0_dp
    else
      v = simpson_rec(f,a,c,fa,flm,fm,left,0.5_dp*tol,depth-1) + &
          simpson_rec(f,c,b,fm,frm,fb,right,0.5_dp*tol,depth-1)
    end if
  end function simpson_rec

  function bessel_k_integral(nu, z) result(k)
    real(dp), intent(in) :: nu, z
    real(dp) :: k, old, upper, h, t, y
    integer :: n, i
    if (z <= 0.0_dp) then
      k = huge_dp
      return
    end if
    upper = max(12.0_dp, log(80.0_dp/z + 1.0_dp) + 4.0_dp)
    n = 128
    old = -1.0_dp
    do
      if (mod(n,2) /= 0) n = n + 1
      h = upper/real(n,dp)
      k = exp(-z)
      do i = 1, n-1
        t = real(i,dp)*h
        y = exp(-z*cosh(t))*cosh(nu*t)
        if (mod(i,2) == 0) then
          k = k + 2.0_dp*y
        else
          k = k + 4.0_dp*y
        end if
      end do
      k = h*k/3.0_dp
      if (old > 0.0_dp) then
        if (abs(k-old) <= 2.0e-11_dp*max(1.0_dp,abs(k))) exit
      end if
      if (n >= 1048576) exit
      old = k
      n = 2*n
    end do
  end function bessel_k_integral

end module actuar_special
