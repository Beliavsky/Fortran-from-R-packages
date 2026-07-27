! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_math
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use peerperformance_kinds, only: dp, pi, sqrt_two
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile, student_t_cdf
  public :: two_sided_normal_pvalue, two_sided_t_pvalue
  public :: beta_inc, finite_value, missing_value, clamp_probability
  public :: set_random_seed, random_integer, random_geometric

contains

  pure elemental real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  pure elemental real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp * erfc(-x / sqrt_two)
  end function normal_cdf

  pure elemental real(dp) function normal_quantile(p) result(x)
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
       4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    real(dp) :: q, r, e

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
    e = normal_cdf(x)-p
    x = x-e/max(normal_pdf(x),tiny(1.0_dp))
  end function normal_quantile

  pure real(dp) function beta_continued_fraction(a, b, x) result(value)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: max_iter = 10000
    real(dp), parameter :: eps = 2.0e-15_dp
    real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
    integer :: m, m2
    real(dp) :: aa, c, d, delta, h, qab, qam, qap
    qab = a+b
    qap = a+1.0_dp
    qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, max_iter
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
      delta = d*c
      h = h*delta
      if (abs(delta-1.0_dp) <= eps) exit
    end do
    value = h
  end function beta_continued_fraction

  pure real(dp) function beta_inc(a, b, x) result(value)
    real(dp), intent(in) :: a, b, x
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      value = 0.0_dp
    else if (x <= 0.0_dp) then
      value = 0.0_dp
    else if (x >= 1.0_dp) then
      value = 1.0_dp
    else
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        value = bt*beta_continued_fraction(a,b,x)/a
      else
        value = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
      end if
      value = clamp_probability(value)
    end if
  end function beta_inc

  pure elemental real(dp) function student_t_cdf(x, df) result(value)
    real(dp), intent(in) :: x, df
    real(dp) :: z, ib
    if (df <= 0.0_dp) then
      value = missing_value()
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      value = 0.5_dp
      return
    end if
    z = df/(df+x*x)
    ib = beta_inc(0.5_dp*df,0.5_dp,z)
    if (x > 0.0_dp) then
      value = 1.0_dp-0.5_dp*ib
    else
      value = 0.5_dp*ib
    end if
    value = clamp_probability(value)
  end function student_t_cdf

  pure elemental real(dp) function two_sided_normal_pvalue(tstat) result(value)
    real(dp), intent(in) :: tstat
    if (.not. ieee_is_finite(tstat)) then
      value = missing_value()
    else
      value = clamp_probability(2.0_dp*normal_cdf(-abs(tstat)))
    end if
  end function two_sided_normal_pvalue

  pure elemental real(dp) function two_sided_t_pvalue(tstat, df) result(value)
    real(dp), intent(in) :: tstat, df
    if (.not. ieee_is_finite(tstat) .or. df <= 0.0_dp) then
      value = missing_value()
    else
      value = clamp_probability(2.0_dp*(1.0_dp-student_t_cdf(abs(tstat),df)))
    end if
  end function two_sided_t_pvalue

  pure elemental logical function finite_value(x) result(ok)
    real(dp), intent(in) :: x
    ok = ieee_is_finite(x)
  end function finite_value

  pure elemental real(dp) function missing_value() result(value)
    value = ieee_value(0.0_dp,ieee_quiet_nan)
  end function missing_value

  pure elemental real(dp) function clamp_probability(x) result(value)
    real(dp), intent(in) :: x
    value = min(1.0_dp,max(0.0_dp,x))
  end function clamp_probability

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: values(:)
    call random_seed(size=n)
    allocate(values(n))
    do i = 1, n
      values(i) = modulo(seed+104729*i,huge(1)-1)
      if (values(i) <= 0) values(i) = i
    end do
    call random_seed(put=values)
  end subroutine set_random_seed

  integer function random_integer(n) result(value)
    integer, intent(in) :: n
    real(dp) :: u
    if (n <= 1) then
      value = 1
      return
    end if
    call random_number(u)
    value = min(n,1+int(u*real(n,dp)))
  end function random_integer

  integer function random_geometric(probability) result(value)
    real(dp), intent(in) :: probability
    real(dp) :: u
    if (probability >= 1.0_dp) then
      value = 0
      return
    end if
    call random_number(u)
    u = max(u,tiny(1.0_dp))
    value = int(log(u)/log(1.0_dp-probability))
  end function random_geometric

end module peerperformance_math
