! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use fextremes_kinds, only: dp, pi
  implicit none
  private
  public :: sort_real, order_real, mean_value, variance_value, quantile_type1
  public :: normal_cdf, normal_quantile, chi_square1_quantile, autocorrelation
  public :: linear_fit, finite_mean_sd, nan_value, clamp
contains
  pure real(dp) function nan_value() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_value

  pure real(dp) function clamp(x, lo, hi) result(y)
    real(dp), intent(in) :: x, lo, hi
    y = min(max(x, lo), hi)
  end function clamp

  subroutine sort_real(x, ascending)
    real(dp), intent(inout) :: x(:)
    logical, intent(in), optional :: ascending
    integer :: i, j
    real(dp) :: key
    logical :: asc
    asc = .true.
    if (present(ascending)) asc = ascending
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      if (asc) then
        do while (j >= 1)
          if (x(j) <= key) exit
          x(j + 1) = x(j)
          j = j - 1
        end do
      else
        do while (j >= 1)
          if (x(j) >= key) exit
          x(j + 1) = x(j)
          j = j - 1
        end do
      end if
      x(j + 1) = key
    end do
  end subroutine sort_real

  subroutine order_real(x, idx, ascending)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: idx(size(x))
    logical, intent(in), optional :: ascending
    integer :: i, j, key
    logical :: asc
    asc = .true.
    if (present(ascending)) asc = ascending
    do i=1,size(x)
      idx(i)=i
    end do
    do i = 2, size(x)
      key = idx(i)
      j = i - 1
      if (asc) then
        do while (j >= 1)
          if (x(idx(j)) <= x(key)) exit
          idx(j + 1) = idx(j)
          j = j - 1
        end do
      else
        do while (j >= 1)
          if (x(idx(j)) >= x(key)) exit
          idx(j + 1) = idx(j)
          j = j - 1
        end do
      end if
      idx(j + 1) = key
    end do
  end subroutine order_real

  pure real(dp) function mean_value(x) result(m)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      m = nan_value()
    else
      m = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure real(dp) function variance_value(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      v = nan_value()
    else
      m = mean_value(x)
      v = sum((x - m)**2) / real(size(x) - 1, dp)
    end if
  end function variance_value

  real(dp) function quantile_type1(x, p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp), allocatable :: work(:)
    integer :: k
    if (size(x) == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      q = nan_value()
      return
    end if
    allocate(work(size(x)))
    work = x
    call sort_real(work)
    if (p <= 0.0_dp) then
      q = work(1)
    else
      k = max(1, min(size(work), ceiling(p * real(size(work), dp))))
      q = work(k)
    end if
  end function quantile_type1

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
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

  pure real(dp) function chi_square1_quantile(p) result(q)
    real(dp), intent(in) :: p
    q = normal_quantile(0.5_dp * (1.0_dp + p))**2
  end function chi_square1_quantile

  subroutine autocorrelation(x, max_lag, acf)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: max_lag
    real(dp), intent(out) :: acf(0:max_lag)
    real(dp) :: m, denom
    integer :: k, n
    n = size(x)
    m = mean_value(x)
    denom = sum((x - m)**2)
    acf = nan_value()
    if (denom <= 0.0_dp) return
    acf(0) = 1.0_dp
    do k = 1, min(max_lag, n - 1)
      acf(k) = sum((x(1:n-k)-m)*(x(1+k:n)-m)) / denom
    end do
  end subroutine autocorrelation

  subroutine linear_fit(x, y, intercept, slope, ok)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: intercept, slope
    logical, intent(out) :: ok
    real(dp) :: mx, my, den
    if (size(x) /= size(y) .or. size(x) < 2) then
      ok = .false.; intercept = nan_value(); slope = nan_value(); return
    end if
    mx = mean_value(x); my = mean_value(y)
    den = sum((x-mx)**2)
    if (den <= epsilon(1.0_dp)) then
      ok = .false.; intercept = nan_value(); slope = nan_value(); return
    end if
    slope = sum((x-mx)*(y-my)) / den
    intercept = my - slope*mx
    ok = .true.
  end subroutine linear_fit

  subroutine finite_mean_sd(x, m, sd, nvalid)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: m, sd
    integer, intent(out) :: nvalid
    integer :: i
    real(dp) :: s, ss
    nvalid = 0; s = 0.0_dp; ss = 0.0_dp
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        nvalid = nvalid + 1
        s = s + x(i); ss = ss + x(i)*x(i)
      end if
    end do
    if (nvalid == 0) then
      m = nan_value(); sd = nan_value()
    else if (nvalid == 1) then
      m = s; sd = nan_value()
    else
      m = s/real(nvalid,dp)
      sd = sqrt(max(0.0_dp,(ss-real(nvalid,dp)*m*m)/real(nvalid-1,dp)))
    end if
  end subroutine finite_mean_sd
end module fextremes_stats
