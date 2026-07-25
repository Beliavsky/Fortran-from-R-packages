! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae_special
  use fcopulae_kinds, only : dp, pi, sqrt_two
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile, student_t_pdf, student_t_cdf
  public :: student_t_quantile, regularized_beta
contains
  elemental function normal_pdf(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: f
    f = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp*erfc(-x/sqrt_two)
  end function normal_cdf

  pure function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  function beta_cf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    real(dp) :: cf
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    real(dp), parameter :: fpmin = 1.0e-300_dp, eps = 3.0e-14_dp
    qab = a+b; qap = a+1.0_dp; qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, 300
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d; if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c; if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d; h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp + aa*d; if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c; if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d; del = d*c; h = h*del
      if (abs(del-1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_cf

  function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: p, bt
    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (x >= 1.0_dp) then
      p = 1.0_dp
    else
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        p = bt*beta_cf(a,b,x)/a
      else
        p = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
      end if
    end if
  end function regularized_beta

  elemental function student_t_pdf(x, nu) result(f)
    real(dp), intent(in) :: x, nu
    real(dp) :: f
    if (nu > 1.0e12_dp) then
      f = normal_pdf(x)
    else
      f = exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)) / &
          sqrt(nu*pi) * (1.0_dp+x*x/nu)**(-0.5_dp*(nu+1.0_dp))
    end if
  end function student_t_pdf

  function student_t_cdf(x, nu) result(p)
    real(dp), intent(in) :: x, nu
    real(dp) :: p, z
    if (nu > 1.0e12_dp) then
      p = normal_cdf(x)
    else if (abs(x) <= tiny(1.0_dp)) then
      p = 0.5_dp
    else
      z = nu/(nu+x*x)
      p = 0.5_dp*regularized_beta(z, 0.5_dp*nu, 0.5_dp)
      if (x > 0.0_dp) p = 1.0_dp-p
    end if
  end function student_t_cdf

  function student_t_quantile(p, nu) result(x)
    real(dp), intent(in) :: p, nu
    real(dp) :: x, lo, hi, mid
    integer :: i
    if (nu > 1.0e12_dp) then
      x = normal_quantile(p)
      return
    end if
    if (p <= 0.0_dp) then; x = -huge(1.0_dp); return; end if
    if (p >= 1.0_dp) then; x = huge(1.0_dp); return; end if
    lo = -1.0_dp; hi = 1.0_dp
    do while (student_t_cdf(lo,nu) > p); lo = 2.0_dp*lo; end do
    do while (student_t_cdf(hi,nu) < p); hi = 2.0_dp*hi; end do
    do i = 1, 100
      mid = 0.5_dp*(lo+hi)
      if (student_t_cdf(mid,nu) < p) then; lo = mid; else; hi = mid; end if
    end do
    x = 0.5_dp*(lo+hi)
  end function student_t_quantile
end module fcopulae_special
