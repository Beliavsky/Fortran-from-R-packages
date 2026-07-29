! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_special
  use copula_kinds, only : dp, pi, sqrt_two
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: student_pdf, student_cdf, student_quantile
  public :: regularized_beta, regularized_gamma_p, gamma_quantile
  public :: log1mexp, log1pexp, log_one_plus, exp_minus_one, debye1
  public :: bivariate_normal_cdf, bivariate_student_cdf
contains
  elemental real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  elemental real(dp) function normal_cdf(x) result(v)
    real(dp), intent(in) :: x
    v = 0.5_dp*erfc(-x/sqrt_two)
  end function normal_cdf

  elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [-3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, &
      2.506628277459239_dp]
    real(dp), parameter :: b(5) = [-5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp]
    real(dp), parameter :: c(6) = [-7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, &
      2.938163982698783_dp]
    real(dp), parameter :: d(4) = [7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp]
    real(dp) :: q, r, e, u
    integer :: i
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= 0.97575_dp) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    do i = 1, 2
      e = normal_cdf(x)-p
      u = e/max(normal_pdf(x), tiny(1.0_dp))
      x = x-u/(1.0_dp+0.5_dp*x*u)
    end do
  end function normal_quantile

  pure real(dp) function beta_cf(a, b, x) result(h)
    real(dp), intent(in) :: a, b, x
    real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
    real(dp) :: qab, qap, qam, c, d, aa, del
    integer :: m, m2
    qab = a+b
    qap = a+1.0_dp
    qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, 400
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) <= eps) exit
    end do
  end function beta_cf

  pure real(dp) function regularized_beta(x, a, b) result(v)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. x <= 0.0_dp) then
      v = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      v = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      v = bt*beta_cf(a,b,x)/a
    else
      v = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    v = min(1.0_dp,max(0.0_dp,v))
  end function regularized_beta

  pure real(dp) function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
    real(dp) :: sumv, del, ap, b, c, d, h, an
    integer :: n
    if (a <= 0.0_dp .or. x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x < a+1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, 600
        ap = ap+1.0_dp
        del = del*x/ap
        sumv = sumv+del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x+1.0_dp-a
      c = 1.0_dp/tiny(1.0_dp)
      d = 1.0_dp/b
      h = d
      do n = 1, 600
        an = -real(n,dp)*(real(n,dp)-a)
        b = b+2.0_dp
        d = an*d+b
        if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
        c = b+an/c
        if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      p = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p = min(1.0_dp,max(0.0_dp,p))
  end function regularized_gamma_p

  real(dp) function gamma_quantile(p, shape, scale) result(x)
    real(dp), intent(in) :: p, shape, scale
    real(dp) :: lo, hi, mid
    integer :: i
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = 0.0_dp
    hi = max(scale,shape*scale)
    do while (regularized_gamma_p(shape,hi/scale) < p)
      hi = 2.0_dp*hi
    end do
    do i = 1, 100
      mid = 0.5_dp*(lo+hi)
      if (regularized_gamma_p(shape,mid/scale) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function gamma_quantile

  elemental real(dp) function student_pdf(x, df) result(v)
    real(dp), intent(in) :: x, df
    if (df > 1.0e12_dp) then
      v = normal_pdf(x)
    else
      v = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df)- &
        0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log(1.0_dp+x*x/df))
    end if
  end function student_pdf

  elemental real(dp) function student_cdf(x, df) result(v)
    real(dp), intent(in) :: x, df
    real(dp) :: z, b
    if (df > 1.0e12_dp) then
      v = normal_cdf(x)
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      v = 0.5_dp
      return
    end if
    z = df/(df+x*x)
    b = regularized_beta(z,0.5_dp*df,0.5_dp)
    if (x > 0.0_dp) then
      v = 1.0_dp-0.5_dp*b
    else
      v = 0.5_dp*b
    end if
  end function student_cdf

  real(dp) function student_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: lo, hi, mid
    integer :: i
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = -1.0_dp
    hi = 1.0_dp
    do while (student_cdf(lo,df) > p)
      lo = 2.0_dp*lo
    end do
    do while (student_cdf(hi,df) < p)
      hi = 2.0_dp*hi
    end do
    do i = 1, 120
      mid = 0.5_dp*(lo+hi)
      if (student_cdf(mid,df) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function student_quantile


  elemental real(dp) function exp_minus_one(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-5_dp) then
      v = x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x*(1.0_dp/24.0_dp+x/120.0_dp))))
    else
      v = exp(x)-1.0_dp
    end if
  end function exp_minus_one

  elemental real(dp) function log_one_plus(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: term, total
    integer :: k
    if (abs(x) < 1.0e-4_dp) then
      term = x
      total = 0.0_dp
      do k = 1, 12
        if (mod(k,2) == 1) then
          total = total+term/real(k,dp)
        else
          total = total-term/real(k,dp)
        end if
        term = term*x
      end do
      v = total
    else
      v = log(1.0_dp+x)
    end if
  end function log_one_plus

  elemental real(dp) function log1mexp(x) result(v)
    real(dp), intent(in) :: x
    if (x > -log(2.0_dp)) then
      v = log(-exp_minus_one(x))
    else
      v = log_one_plus(-exp(x))
    end if
  end function log1mexp

  elemental real(dp) function log1pexp(x) result(v)
    real(dp), intent(in) :: x
    if (x > 33.3_dp) then
      v = x
    else if (x < -37.0_dp) then
      v = exp(x)
    else
      v = log_one_plus(exp(x))
    end if
  end function log1pexp


  real(dp) function debye1(theta) result(value)
    real(dp), intent(in) :: theta
    integer, parameter :: n = 4096
    real(dp) :: h, x, total, fx
    integer :: i
    if (abs(theta) < 1.0e-5_dp) then
      value = 1.0_dp-theta/4.0_dp+theta*theta/36.0_dp-theta**4/3600.0_dp
      return
    end if
    h = theta/real(n,dp)
    total = 1.0_dp+theta/exp_minus_one(theta)
    do i = 1, n-1
      x = real(i,dp)*h
      if (abs(x) < 1.0e-7_dp) then
        fx = 1.0_dp-x/2.0_dp+x*x/12.0_dp
      else
        fx = x/exp_minus_one(x)
      end if
      if (mod(i,2) == 0) then
        total = total+2.0_dp*fx
      else
        total = total+4.0_dp*fx
      end if
    end do
    value = h*total/(3.0_dp*theta)
  end function debye1

  real(dp) function bivariate_normal_cdf(x, y, rho) result(value)
    real(dp), intent(in) :: x, y, rho
    integer, parameter :: n = 4096
    real(dp) :: lower, upper, h, z, total, fz
    integer :: i
    if (rho >= 1.0_dp-1.0e-14_dp) then
      value = normal_cdf(min(x,y))
      return
    else if (rho <= -1.0_dp+1.0e-14_dp) then
      value = max(0.0_dp,normal_cdf(x)-normal_cdf(-y))
      return
    end if
    if (x < -10.0_dp .or. y < -10.0_dp) then
      value = 0.0_dp
      return
    else if (x > 10.0_dp) then
      value = normal_cdf(y)
      return
    else if (y > 10.0_dp) then
      value = normal_cdf(x)
      return
    end if
    lower = -10.0_dp
    upper = x
    h = (upper-lower)/real(n,dp)
    total = normal_pdf(lower)*normal_cdf((y-rho*lower)/sqrt(1.0_dp-rho*rho)) + &
      normal_pdf(upper)*normal_cdf((y-rho*upper)/sqrt(1.0_dp-rho*rho))
    do i = 1, n-1
      z = lower+real(i,dp)*h
      fz = normal_pdf(z)*normal_cdf((y-rho*z)/sqrt(1.0_dp-rho*rho))
      if (mod(i,2) == 0) then
        total = total+2.0_dp*fz
      else
        total = total+4.0_dp*fz
      end if
    end do
    value = h*total/3.0_dp
    value = min(1.0_dp,max(0.0_dp,value))
  end function bivariate_normal_cdf

  real(dp) function bivariate_student_cdf(x, y, rho, df) result(value)
    real(dp), intent(in) :: x, y, rho, df
    integer, parameter :: n = 8192
    real(dp) :: lower, upper, h, z, total, fz, scale, arg
    integer :: i
    if (rho >= 1.0_dp-1.0e-14_dp) then
      value = student_cdf(min(x,y),df)
      return
    else if (rho <= -1.0_dp+1.0e-14_dp) then
      value = max(0.0_dp,student_cdf(x,df)-student_cdf(-y,df))
      return
    end if
    lower = -max(30.0_dp,10.0_dp*sqrt(df/max(df-2.0_dp,0.2_dp)))
    upper = x
    if (upper <= lower) then
      value = 0.0_dp
      return
    end if
    h = (upper-lower)/real(n,dp)
    scale = sqrt((df+lower*lower)*(1.0_dp-rho*rho)/(df+1.0_dp))
    arg = (y-rho*lower)/scale
    total = student_pdf(lower,df)*student_cdf(arg,df+1.0_dp)
    scale = sqrt((df+upper*upper)*(1.0_dp-rho*rho)/(df+1.0_dp))
    arg = (y-rho*upper)/scale
    total = total+student_pdf(upper,df)*student_cdf(arg,df+1.0_dp)
    do i = 1, n-1
      z = lower+real(i,dp)*h
      scale = sqrt((df+z*z)*(1.0_dp-rho*rho)/(df+1.0_dp))
      arg = (y-rho*z)/scale
      fz = student_pdf(z,df)*student_cdf(arg,df+1.0_dp)
      if (mod(i,2) == 0) then
        total = total+2.0_dp*fz
      else
        total = total+4.0_dp*fz
      end if
    end do
    value = h*total/3.0_dp
    value = min(1.0_dp,max(0.0_dp,value))
  end function bivariate_student_cdf
end module copula_special
