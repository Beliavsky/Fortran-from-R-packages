! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_special
  use mvtnorm_kinds, only : dp, pi, sqrt_two
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: log_gamma_dp, regularized_beta, student_t_cdf, student_t_quantile
  public :: regularized_gamma_p, chi_square_cdf, chi_square_quantile
  public :: log1mexp, logsumexp2

contains

  elemental real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  elemental real(dp) function normal_cdf(x) result(v)
    real(dp), intent(in) :: x
    v = 0.5_dp*erfc(-x/sqrt_two)
  end function normal_cdf

  elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e+01_dp,  2.209460984245205e+02_dp, &
      -2.759285104469687e+02_dp,  1.383577518672690e+02_dp, &
      -3.066479806614716e+01_dp,  2.506628277459239e+00_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e+01_dp,  1.615858368580409e+02_dp, &
      -1.556989798598866e+02_dp,  6.680131188771972e+01_dp, &
      -1.328068155288572e+01_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
      -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
       4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp,  3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp,  3.754408661907416e+00_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    real(dp) :: q, r, e, u
    integer :: i

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if

    if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if

    ! Two Halley refinements are enough to recover full double precision.
    do i = 1, 2
      e = normal_cdf(x)-p
      u = e/max(normal_pdf(x), tiny(1.0_dp))
      x = x-u/(1.0_dp+0.5_dp*x*u)
    end do
  end function normal_quantile

  elemental real(dp) function log_gamma_dp(x) result(v)
    real(dp), intent(in) :: x
    v = log_gamma(x)
  end function log_gamma_dp

  real(dp) function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 300
    real(dp), parameter :: eps = 4.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
    real(dp) :: qab, qap, qam, c, d, h, aa, del
    integer :: m, m2

    qab = a+b
    qap = a+1.0_dp
    qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, maxit
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
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
  end function beta_continued_fraction

  real(dp) function regularized_beta(x, a, b) result(v)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      v = 0.0_dp
      return
    end if
    if (x <= 0.0_dp) then
      v = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      v = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log1p_local(-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      v = bt*beta_continued_fraction(a,b,x)/a
    else
      v = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
    end if
    v = min(1.0_dp,max(0.0_dp,v))
  end function regularized_beta

  real(dp) function student_t_cdf(x, df) result(v)
    real(dp), intent(in) :: x, df
    real(dp) :: z, ib
    if (df <= 0.0_dp) then
      v = normal_cdf(x)
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      v = 0.5_dp
      return
    end if
    z = df/(df+x*x)
    ib = regularized_beta(z,0.5_dp*df,0.5_dp)
    if (x > 0.0_dp) then
      v = 1.0_dp-0.5_dp*ib
    else
      v = 0.5_dp*ib
    end if
  end function student_t_cdf

  real(dp) function student_t_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: lo, hi, mid, f, pdf
    integer :: i
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (df <= 0.0_dp) then
      x = normal_quantile(p)
      return
    end if
    x = normal_quantile(p)*sqrt(max(df/(df-2.0_dp),1.0_dp))
    lo = -1.0_dp
    hi = 1.0_dp
    do while (student_t_cdf(lo,df) > p)
      lo = 2.0_dp*lo
    end do
    do while (student_t_cdf(hi,df) < p)
      hi = 2.0_dp*hi
    end do
    x = min(hi,max(lo,x))
    do i = 1, 80
      f = student_t_cdf(x,df)-p
      pdf = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df) &
            -0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log1p_local(x*x/df))
      if (f > 0.0_dp) hi = x
      if (f < 0.0_dp) lo = x
      if (pdf > tiny(1.0_dp)) then
        mid = x-f/pdf
      else
        mid = 0.5_dp*(lo+hi)
      end if
      if (mid <= lo .or. mid >= hi) mid = 0.5_dp*(lo+hi)
      if (abs(mid-x) <= 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x))) exit
      x = mid
    end do
  end function student_t_quantile

  real(dp) function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: maxit = 500
    real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
    real(dp) :: sumv, del, ap, b, c, d, h, an
    integer :: n
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      p = 0.0_dp
      return
    end if
    if (x < a+1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, maxit
        ap = ap+1.0_dp
        del = del*x/ap
        sumv = sumv+del
        if (abs(del) < abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x+1.0_dp-a
      c = 1.0_dp/tiny(1.0_dp)
      d = 1.0_dp/b
      h = d
      do n = 1, maxit
        an = -real(n,dp)*(real(n,dp)-a)
        b = b+2.0_dp
        d = an*d+b
        if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
        c = b+an/c
        if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) < eps) exit
      end do
      p = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p = min(1.0_dp,max(0.0_dp,p))
  end function regularized_gamma_p

  real(dp) function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = regularized_gamma_p(0.5_dp*df,0.5_dp*x)
    end if
  end function chi_square_cdf

  real(dp) function chi_square_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: lo, hi, mid, f, pdf, z
    integer :: i
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    z = normal_quantile(p)
    x = df*(1.0_dp-2.0_dp/(9.0_dp*df)+z*sqrt(2.0_dp/(9.0_dp*df)))**3
    x = max(x, tiny(1.0_dp))
    lo = 0.0_dp
    hi = max(df+10.0_dp*sqrt(2.0_dp*df),2.0_dp*x)
    do while (chi_square_cdf(hi,df) < p)
      hi = 2.0_dp*hi
    end do
    do i = 1, 80
      f = chi_square_cdf(x,df)-p
      pdf = exp((0.5_dp*df-1.0_dp)*log(x)-0.5_dp*x &
            -0.5_dp*df*log(2.0_dp)-log_gamma(0.5_dp*df))
      if (f > 0.0_dp) hi = x
      if (f < 0.0_dp) lo = x
      if (pdf > tiny(1.0_dp)) then
        mid = x-f/pdf
      else
        mid = 0.5_dp*(lo+hi)
      end if
      if (mid <= lo .or. mid >= hi) mid = 0.5_dp*(lo+hi)
      if (abs(mid-x) <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,x)) exit
      x = mid
    end do
  end function chi_square_quantile

  elemental real(dp) function log1mexp(x) result(v)
    real(dp), intent(in) :: x
    if (x >= 0.0_dp) then
      v = -huge(1.0_dp)
    else if (x < -log(2.0_dp)) then
      v = log1p_local(-exp(x))
    else
      v = log(-expm1_local(x))
    end if
  end function log1mexp

  elemental real(dp) function logsumexp2(a,b) result(v)
    real(dp), intent(in) :: a,b
    real(dp) :: m
    m = max(a,b)
    if (m <= -huge(1.0_dp)/2.0_dp) then
      v = m
    else
      v = m+log(exp(a-m)+exp(b-m))
    end if
  end function logsumexp2

  elemental real(dp) function log1p_local(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-4_dp) then
      v = x*(1.0_dp-x*(0.5_dp-x*(1.0_dp/3.0_dp-x*0.25_dp)))
    else
      v = log(1.0_dp+x)
    end if
  end function log1p_local

  elemental real(dp) function expm1_local(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-5_dp) then
      v = x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x/24.0_dp)))
    else
      v = exp(x)-1.0_dp
    end if
  end function expm1_local

end module mvtnorm_special
