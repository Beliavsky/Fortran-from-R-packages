! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

module betategarch_math
  use betategarch_kinds, only : dp
  implicit none
  private

  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: huge_penalty = -1.0e100_dp

  public :: beta_fn, log_beta_fn, signum, all_finite
  public :: invert_matrix, numerical_hessian

  abstract interface
    function scalar_function(x, context) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: context
      real(dp) :: value
    end function scalar_function
  end interface

contains

  pure function log_beta_fn(a, b) result(value)
    real(dp), intent(in) :: a, b
    real(dp) :: value

    value = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
  end function log_beta_fn

  pure function beta_fn(a, b) result(value)
    real(dp), intent(in) :: a, b
    real(dp) :: value

    value = exp(log_beta_fn(a, b))
  end function beta_fn

  pure elemental function signum(x) result(value)
    real(dp), intent(in) :: x
    integer :: value

    if (x > 0.0_dp) then
      value = 1
    else if (x < 0.0_dp) then
      value = -1
    else
      value = 0
    end if
  end function signum

  pure function all_finite(x) result(ok)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp), intent(in) :: x(:)
    logical :: ok

    ok = all(ieee_is_finite(x))
  end function all_finite

  subroutine invert_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: ainv(:, :)
    integer, intent(out) :: info

    real(dp), allocatable :: aug(:, :), row_tmp(:)
    real(dp) :: pivot_abs, factor
    integer :: n, i, j, pivot_row

    n = size(a, 1)
    info = 0
    ainv = 0.0_dp
    if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
      info = -1
      return
    end if

    allocate(aug(n, 2*n), row_tmp(2*n))
    aug(:, 1:n) = a
    aug(:, n+1:2*n) = 0.0_dp
    do i = 1, n
      aug(i, n+i) = 1.0_dp
    end do

    do i = 1, n
      pivot_row = i
      pivot_abs = abs(aug(i, i))
      do j = i + 1, n
        if (abs(aug(j, i)) > pivot_abs) then
          pivot_abs = abs(aug(j, i))
          pivot_row = j
        end if
      end do
      if (pivot_abs <= sqrt(tiny(1.0_dp))) then
        info = i
        return
      end if
      if (pivot_row /= i) then
        row_tmp = aug(i, :)
        aug(i, :) = aug(pivot_row, :)
        aug(pivot_row, :) = row_tmp
      end if

      aug(i, :) = aug(i, :) / aug(i, i)
      do j = 1, n
        if (j /= i) then
          factor = aug(j, i)
          aug(j, :) = aug(j, :) - factor * aug(i, :)
        end if
      end do
    end do

    ainv = aug(:, n+1:2*n)
  end subroutine invert_matrix

  subroutine numerical_hessian(fun, x, context, hessian, info, relative_step)
    procedure(scalar_function) :: fun
    real(dp), intent(in) :: x(:)
    class(*), intent(in) :: context
    real(dp), intent(out) :: hessian(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: relative_step

    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), step(:)
    real(dp) :: f0, hi, hj, rel
    integer :: n, i, j

    n = size(x)
    info = 0
    if (size(hessian, 1) /= n .or. size(hessian, 2) /= n) then
      info = -1
      return
    end if

    rel = epsilon(1.0_dp)**0.25_dp
    if (present(relative_step)) rel = relative_step
    allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n), step(n))
    step = rel * max(1.0_dp, abs(x))
    f0 = fun(x, context)

    do i = 1, n
      hi = step(i)
      xp = x
      xm = x
      xp(i) = xp(i) + hi
      xm(i) = xm(i) - hi
      hessian(i, i) = (fun(xp, context) - 2.0_dp*f0 + fun(xm, context)) / (hi*hi)
      do j = i + 1, n
        hj = step(j)
        xpp = x
        xpm = x
        xmp = x
        xmm = x
        xpp(i) = xpp(i) + hi
        xpp(j) = xpp(j) + hj
        xpm(i) = xpm(i) + hi
        xpm(j) = xpm(j) - hj
        xmp(i) = xmp(i) - hi
        xmp(j) = xmp(j) + hj
        xmm(i) = xmm(i) - hi
        xmm(j) = xmm(j) - hj
        hessian(i, j) = (fun(xpp, context) - fun(xpm, context) - fun(xmp, context) + fun(xmm, context)) / (4.0_dp*hi*hj)
        hessian(j, i) = hessian(i, j)
      end do
    end do
  end subroutine numerical_hessian

end module betategarch_math
