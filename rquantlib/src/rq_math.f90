! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_math
  use rq_kinds, only: dp, pi, sqrt_two_pi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_inv, clamp, nearly_equal
  public :: seed_random, random_uniform, random_normal
  public :: bisect_root, sort_real

  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
  end interface
contains
  pure elemental function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = exp(-0.5_dp*x*x)/sqrt_two_pi
  end function normal_pdf

  pure elemental function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure elemental function clamp(x, lo, hi) result(y)
    real(dp), intent(in) :: x, lo, hi
    real(dp) :: y
    y = min(max(x, lo), hi)
  end function clamp

  pure function nearly_equal(a, b, atol, rtol) result(ok)
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: atol, rtol
    logical :: ok
    real(dp) :: aa, rr
    aa = 1.0e-12_dp
    rr = 1.0e-10_dp
    if (present(atol)) aa = atol
    if (present(rtol)) rr = rtol
    ok = abs(a-b) <= aa + rr*max(abs(a),abs(b))
  end function nearly_equal

  pure elemental function normal_inv(p) result(x)
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
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > phigh) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_inv

  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
      if (put(i) <= 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine seed_random

  function random_uniform() result(u)
    real(dp) :: u
    call random_number(u)
    u = max(u, tiny(1.0_dp))
  end function random_uniform

  function random_normal() result(z)
    real(dp) :: z
    z = sqrt(-2.0_dp*log(random_uniform()))*cos(2.0_dp*pi*random_uniform())
  end function random_normal

  subroutine bisect_root(fun, lo, hi, root, status, tol, max_iter)
    procedure(scalar_function) :: fun
    real(dp), intent(in) :: lo, hi
    real(dp), intent(out) :: root
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    real(dp) :: a, b, c, fa, fb, fc, eps
    integer :: i, nmax
    eps = 1.0e-10_dp
    if (present(tol)) eps = tol
    nmax = 200
    if (present(max_iter)) nmax = max_iter
    a = lo; b = hi; fa = fun(a); fb = fun(b)
    if (ieee_is_nan(fa) .or. ieee_is_nan(fb) .or. fa*fb > 0.0_dp) then
      root = 0.5_dp*(a+b); status = 1; return
    end if
    do i=1,nmax
      c = 0.5_dp*(a+b); fc = fun(c)
      if (abs(fc) <= eps .or. 0.5_dp*(b-a) <= eps*(1.0_dp+abs(c))) then
        root = c; status = 0; return
      end if
      if (fa*fc <= 0.0_dp) then
        b = c; fb = fc
      else
        a = c; fa = fc
      end if
    end do
    root = 0.5_dp*(a+b); status = 2
  end subroutine bisect_root

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while (j>=1)
        if (x(j) <= key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real
end module rq_math
