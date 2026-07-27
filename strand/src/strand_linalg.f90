! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand_linalg
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use strand_kinds, only : dp
  implicit none
  private
  public :: least_squares, weighted_least_squares, correlation

contains

  subroutine least_squares(x, y, beta, ok)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), intent(out) :: beta(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: w(:)

    allocate(w(size(y)))
    w = 1.0_dp
    call weighted_least_squares(x, y, w, beta, ok)
  end subroutine least_squares

  subroutine weighted_least_squares(x, y, weights, beta, ok)
    real(dp), intent(in) :: x(:, :), y(:), weights(:)
    real(dp), intent(out) :: beta(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: a(:, :), yy(:), column_norms(:)
    real(dp), allocatable :: q(:, :), r(:, :), rhs(:), scales(:), solution(:), v(:)
    real(dp) :: norm_v, temp, tolerance, sw
    integer, allocatable :: permutation(:)
    integer :: i, j, k, n, p, pivot_col, rank, temp_index

    n = size(x, 1)
    p = size(x, 2)
    if (size(y) /= n .or. size(weights) /= n .or. size(beta) /= p) then
      error stop 'weighted_least_squares: inconsistent dimensions'
    end if
    beta = 0.0_dp
    if (n == 0 .or. p == 0 .or. any(weights < 0.0_dp)) then
      ok = .false.
      return
    end if

    allocate(a(n, p), yy(n), q(n, p), r(p, p), rhs(p), scales(p), solution(p), v(n))
    allocate(column_norms(p), permutation(p))
    do i = 1, n
      sw = sqrt(weights(i))
      a(i, :) = sw * x(i, :)
      yy(i) = sw * y(i)
    end do
    q = 0.0_dp
    r = 0.0_dp
    rhs = 0.0_dp
    solution = 0.0_dp
    do j = 1, p
      scales(j) = sqrt(sum(a(:, j)**2))
      if (scales(j) > tiny(1.0_dp)) then
        a(:, j) = a(:, j) / scales(j)
      else
        a(:, j) = 0.0_dp
        scales(j) = 1.0_dp
      end if
      column_norms(j) = sum(a(:, j)**2)
      permutation(j) = j
    end do

    tolerance = 100.0_dp * epsilon(1.0_dp) * real(max(n, p), dp)
    rank = 0
    do k = 1, p
      pivot_col = k - 1 + maxloc(column_norms(k:p), dim=1)
      if (pivot_col /= k) then
        call swap_columns(a, k, pivot_col)
        temp = column_norms(k)
        column_norms(k) = column_norms(pivot_col)
        column_norms(pivot_col) = temp
        temp_index = permutation(k)
        permutation(k) = permutation(pivot_col)
        permutation(pivot_col) = temp_index
        do i = 1, k - 1
          temp = r(i, k)
          r(i, k) = r(i, pivot_col)
          r(i, pivot_col) = temp
        end do
      end if

      v = a(:, k)
      norm_v = sqrt(sum(v**2))
      if (norm_v <= tolerance) exit
      rank = k
      q(:, k) = v / norm_v
      r(k, k) = norm_v
      do j = k + 1, p
        r(k, j) = dot_product(q(:, k), a(:, j))
        a(:, j) = a(:, j) - r(k, j) * q(:, k)
        column_norms(j) = sum(a(:, j)**2)
      end do
    end do

    if (rank == 0) then
      ok = .false.
      return
    end if

    rhs(1:rank) = matmul(transpose(q(:, 1:rank)), yy)
    solution(rank) = rhs(rank) / r(rank, rank)
    do i = rank - 1, 1, -1
      solution(i) = (rhs(i) - dot_product(r(i, i + 1:rank), solution(i + 1:rank))) / r(i, i)
    end do
    do i = 1, rank
      beta(permutation(i)) = solution(i) / scales(permutation(i))
    end do
    ok = all(ieee_is_finite(beta))
    if (.not. ok) beta = 0.0_dp
  end subroutine weighted_least_squares

  subroutine swap_columns(a, first, second)
    real(dp), intent(inout) :: a(:, :)
    integer, intent(in) :: first, second
    real(dp), allocatable :: temporary(:)

    allocate(temporary(size(a, 1)))
    temporary = a(:, first)
    a(:, first) = a(:, second)
    a(:, second) = temporary
  end subroutine swap_columns

  function correlation(x, y) result(rho)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: rho
    real(dp) :: xm, ym, sx, sy

    if (size(x) /= size(y) .or. size(x) < 2) then
      rho = 0.0_dp
      return
    end if
    xm = sum(x) / real(size(x), dp)
    ym = sum(y) / real(size(y), dp)
    sx = sqrt(sum((x - xm)**2))
    sy = sqrt(sum((y - ym)**2))
    if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
      rho = 0.0_dp
    else
      rho = dot_product(x - xm, y - ym) / (sx * sy)
    end if
  end function correlation

end module strand_linalg
