! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim_linalg
  use portfoliooptim_kinds, only : dp
  implicit none
  private
  public :: diagonal_matrix, solve_linear, invert_matrix, max_abs

contains

  pure function diagonal_matrix(x) result(a)
    real(dp), intent(in) :: x(:)
    real(dp) :: a(size(x), size(x))
    integer :: i

    a = 0.0_dp
    do i = 1, size(x)
      a(i, i) = x(i)
    end do
  end function diagonal_matrix

  pure real(dp) function max_abs(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = maxval(abs(x))
    end if
  end function max_abs

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aa(:, :), bb(:), row_tmp(:)
    real(dp) :: factor, pivot_tmp, scale
    integer :: n, i, j, p

    n = size(b)
    ok = .false.
    x = 0.0_dp
    if (size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) return
    if (n == 0) then
      ok = .true.
      return
    end if

    allocate(aa(n, n), bb(n), row_tmp(n))
    aa = a
    bb = b
    scale = max(1.0_dp, maxval(abs(aa)))

    do i = 1, n
      p = i - 1 + maxloc(abs(aa(i:n, i)), dim=1)
      if (abs(aa(p, i)) <= 100.0_dp * epsilon(1.0_dp) * scale) return
      if (p /= i) then
        row_tmp = aa(i, :)
        aa(i, :) = aa(p, :)
        aa(p, :) = row_tmp
        pivot_tmp = bb(i)
        bb(i) = bb(p)
        bb(p) = pivot_tmp
      end if
      do j = i + 1, n
        factor = aa(j, i) / aa(i, i)
        aa(j, i) = 0.0_dp
        aa(j, i + 1:n) = aa(j, i + 1:n) - factor * aa(i, i + 1:n)
        bb(j) = bb(j) - factor * bb(i)
      end do
    end do

    do i = n, 1, -1
      if (i < n) then
        x(i) = (bb(i) - dot_product(aa(i, i + 1:n), x(i + 1:n))) / aa(i, i)
      else
        x(i) = bb(i) / aa(i, i)
      end if
    end do
    ok = all(ieee_is_finite_vector(x))
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: ainv(:, :)
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:), col(:)
    logical :: col_ok
    integer :: n, j

    n = size(a, 1)
    ok = .false.
    ainv = 0.0_dp
    if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) return
    allocate(e(n), col(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, col, col_ok)
      if (.not. col_ok) return
      ainv(:, j) = col
    end do
    ok = .true.
  end subroutine invert_matrix

  pure function ieee_is_finite_vector(x) result(mask)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp), intent(in) :: x(:)
    logical :: mask(size(x))
    integer :: i
    do i = 1, size(x)
      mask(i) = ieee_is_finite(x(i))
    end do
  end function ieee_is_finite_vector

end module portfoliooptim_linalg
