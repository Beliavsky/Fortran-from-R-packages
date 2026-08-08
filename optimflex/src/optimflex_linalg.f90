! SPDX-License-Identifier: MIT
module optimflex_linalg
  use optimflex_types, only : dp
  implicit none
  private
  public :: cholesky_lower, is_pd_fast, solve_linear, solve_spd
  public :: inverse_matrix, symmetrize, vecnorm2, maxabs

contains

  pure real(dp) function vecnorm2(x) result(v)
    real(dp), intent(in) :: x(:)
    v = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vecnorm2

  pure real(dp) function maxabs(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = 0.0_dp
    else
      v = maxval(abs(x))
    end if
  end function maxabs

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    a = 0.5_dp * (a + transpose(a))
  end subroutine symmetrize

  subroutine cholesky_lower(a, l, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: i, j, k, n
    real(dp) :: s, scale, tol

    n = size(a, 1)
    l = 0.0_dp
    ok = .false.
    if (size(a, 2) /= n .or. size(l,1) /= n .or. size(l,2) /= n) return
    scale = max(1.0_dp, maxval(abs(a)))
    tol = 100.0_dp * epsilon(1.0_dp) * scale
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j - 1
          s = s - l(i,k) * l(j,k)
        end do
        if (i == j) then
          if (s <= tol) return
          l(i,j) = sqrt(s)
        else
          if (abs(l(j,j)) <= tiny(1.0_dp)) return
          l(i,j) = s / l(j,j)
        end if
      end do
    end do
    ok = .true.
  end subroutine cholesky_lower

  logical function is_pd_fast(a) result(ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: l(:,:)
    if (size(a,1) /= size(a,2) .or. size(a,1) == 0) then
      ok = .false.
      return
    end if
    allocate(l(size(a,1), size(a,2)))
    call cholesky_lower(a, l, ok)
  end function is_pd_fast

  subroutine solve_spd(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: l(:,:), y(:)
    integer :: i, j, n

    n = size(b)
    ok = .false.
    x = 0.0_dp
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
    allocate(l(n,n), y(n))
    call cholesky_lower(a, l, ok)
    if (.not. ok) return
    do i = 1, n
      y(i) = b(i)
      do j = 1, i - 1
        y(i) = y(i) - l(i,j) * y(j)
      end do
      y(i) = y(i) / l(i,i)
    end do
    do i = n, 1, -1
      x(i) = y(i)
      do j = i + 1, n
        x(i) = x(i) - l(j,i) * x(j)
      end do
      x(i) = x(i) / l(i,i)
    end do
    ok = .true.
  end subroutine solve_spd

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: m(:,:), rhs(:), rowtmp(:)
    real(dp) :: pivot, fac, scale, tol
    integer :: i, j, k, p, n

    n = size(b)
    ok = .false.
    x = 0.0_dp
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
    allocate(m(n,n), rhs(n), rowtmp(n))
    m = a
    rhs = b
    scale = max(1.0_dp, maxval(abs(m)))
    tol = 100.0_dp * epsilon(1.0_dp) * scale
    do k = 1, n - 1
      p = k - 1 + maxloc(abs(m(k:n,k)), dim=1)
      if (abs(m(p,k)) <= tol) return
      if (p /= k) then
        rowtmp = m(k,:)
        m(k,:) = m(p,:)
        m(p,:) = rowtmp
        pivot = rhs(k)
        rhs(k) = rhs(p)
        rhs(p) = pivot
      end if
      do i = k + 1, n
        fac = m(i,k) / m(k,k)
        m(i,k) = 0.0_dp
        do j = k + 1, n
          m(i,j) = m(i,j) - fac * m(k,j)
        end do
        rhs(i) = rhs(i) - fac * rhs(k)
      end do
    end do
    if (abs(m(n,n)) <= tol) return
    do i = n, 1, -1
      pivot = rhs(i)
      do j = i + 1, n
        pivot = pivot - m(i,j) * x(j)
      end do
      if (abs(m(i,i)) <= tol) return
      x(i) = pivot / m(i,i)
    end do
    ok = .true.
  end subroutine solve_linear

  subroutine inverse_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:), col(:)
    integer :: j, n
    logical :: col_ok

    n = size(a,1)
    ok = .false.
    ainv = 0.0_dp
    if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) return
    allocate(e(n), col(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, col, col_ok)
      if (.not. col_ok) return
      ainv(:,j) = col
    end do
    ok = .true.
  end subroutine inverse_matrix

end module optimflex_linalg
