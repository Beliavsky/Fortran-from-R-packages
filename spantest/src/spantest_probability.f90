! SPDX-License-Identifier: GPL-3.0-only
module spantest_probability
  use spantest_kinds, only : dp
  implicit none
  private
  public :: normal_cdf, normal_quantile, student_t_cdf, f_upper_tail

contains

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
      -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
      -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
      -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
      -1.328068155288572e+01_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
      -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
       4.374664141464968e+00_dp, 2.938163982698783e+00_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
    real(dp), parameter :: plow = 0.02425_dp
    real(dp), parameter :: phigh = 1.0_dp - plow
    real(dp) :: q, r

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 300
    real(dp), parameter :: eps = 3.0e-14_dp
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, maxit
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
           ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
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

  pure real(dp) function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
             a*log(x) + b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      p = bt*beta_continued_fraction(a,b,x)/a
    else
      p = 1.0_dp - bt*beta_continued_fraction(b,a,1.0_dp-x)/b
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function regularized_beta

  pure real(dp) function student_t_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    real(dp) :: z, ib
    if (df <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      p = 0.5_dp
      return
    end if
    z = df / (df + x*x)
    ib = regularized_beta(z, 0.5_dp*df, 0.5_dp)
    if (x > 0.0_dp) then
      p = 1.0_dp - 0.5_dp*ib
    else
      p = 0.5_dp*ib
    end if
  end function student_t_cdf

  pure real(dp) function f_upper_tail(x, df1, df2) result(p)
    real(dp), intent(in) :: x, df1, df2
    real(dp) :: z
    if (x < 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
      p = 1.0_dp
      return
    end if
    z = df2 / (df2 + df1*x)
    p = regularized_beta(z, 0.5_dp*df2, 0.5_dp*df1)
  end function f_upper_tail

end module spantest_probability
