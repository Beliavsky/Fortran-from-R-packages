! SPDX-License-Identifier: GPL-2.0-only
module optimx_linalg
  use optimx_kinds, only: dp
  implicit none
  private
  public :: vector_norm, eye, outer, solve_linear, cholesky_pd, symmetric_min_bound

contains
  pure real(dp) function vector_norm(x) result(v)
    real(dp), intent(in) :: x(:)
    v = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vector_norm

  pure function eye(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function eye

  pure function outer(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: i
    do i = 1, size(x)
      a(i, :) = x(i) * y
    end do
  end function outer

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: m(:, :), rhs(:), row(:)
    real(dp) :: pivot, factor
    integer :: n, i, j, k, p
    n = size(b)
    status = 1
    x = 0.0_dp
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
    allocate(m(n,n), rhs(n), row(n))
    m = a
    rhs = b
    do k = 1, n - 1
      p = k
      do i = k + 1, n
        if (abs(m(i,k)) > abs(m(p,k))) p = i
      end do
      if (abs(m(p,k)) <= tiny(1.0_dp)) return
      if (p /= k) then
        row = m(k,:); m(k,:) = m(p,:); m(p,:) = row
        pivot = rhs(k); rhs(k) = rhs(p); rhs(p) = pivot
      end if
      do i = k + 1, n
        factor = m(i,k) / m(k,k)
        m(i,k) = 0.0_dp
        do j = k + 1, n
          m(i,j) = m(i,j) - factor * m(k,j)
        end do
        rhs(i) = rhs(i) - factor * rhs(k)
      end do
    end do
    if (abs(m(n,n)) <= tiny(1.0_dp)) return
    x(n) = rhs(n) / m(n,n)
    do i = n - 1, 1, -1
      x(i) = (rhs(i) - dot_product(m(i,i+1:n), x(i+1:n))) / m(i,i)
    end do
    status = 0
  end subroutine solve_linear

  subroutine cholesky_pd(a, l, status, tol)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: l(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tol
    real(dp) :: s, threshold
    integer :: n, i, j, k
    n = size(a,1)
    l = 0.0_dp
    status = 1
    if (size(a,2) /= n .or. any(shape(l) /= [n,n])) return
    threshold = 1.0e-12_dp
    if (present(tol)) threshold = tol
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j - 1
          s = s - l(i,k) * l(j,k)
        end do
        if (i == j) then
          if (s <= threshold) return
          l(i,j) = sqrt(s)
        else
          l(i,j) = s / l(j,j)
        end if
      end do
    end do
    status = 0
  end subroutine cholesky_pd

  pure real(dp) function symmetric_min_bound(a) result(v)
    real(dp), intent(in) :: a(:, :)
    integer :: i, j, n
    real(dp) :: radius
    n = size(a,1)
    v = huge(1.0_dp)
    if (size(a,2) /= n) then
      v = -huge(1.0_dp)
      return
    end if
    do i = 1, n
      radius = 0.0_dp
      do j = 1, n
        if (j /= i) radius = radius + abs(a(i,j))
      end do
      v = min(v, a(i,i) - radius)
    end do
  end function symmetric_min_bound
end module optimx_linalg
