! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_probability
  use fportfolio_kinds, only: dp, pi
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile
contains
  pure real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: q, r
    real(dp), parameter :: a(6) = [ -3.969683028665376d1, 2.209460984245205d2, &
      -2.759285104469687d2, 1.383577518672690d2, -3.066479806614716d1, 2.506628277459239d0 ]
    real(dp), parameter :: b(5) = [ -5.447609879822406d1, 1.615858368580409d2, &
      -1.556989798598866d2, 6.680131188771972d1, -1.328068155288572d1 ]
    real(dp), parameter :: c(6) = [ -7.784894002430293d-3, -3.223964580411365d-1, &
      -2.400758277161838d0, -2.549732539343734d0, 4.374664141464968d0, 2.938163982698783d0 ]
    real(dp), parameter :: d(4) = [ 7.784695709041462d-3, 3.224671290700398d-1, &
      2.445134137142996d0, 3.754408661907416d0 ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
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
    if (abs(x) < huge(1.0_dp)/2.0_dp) then
      x = x - (normal_cdf(x)-p)/max(normal_pdf(x), tiny(1.0_dp))
    end if
  end function normal_quantile
end module fportfolio_probability
