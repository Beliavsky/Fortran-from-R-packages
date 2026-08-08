! SPDX-License-Identifier: Apache-2.0
module psqn_linalg
  use psqn_types, only : dp
  implicit none
  private
  public :: packed_size, packed_index, sym_packed_matvec, sym_packed_matvec_indexed
  public :: rank_one_update_packed, bfgs_inverse_update_packed
  public :: kahan_add, packed_to_dense, dense_to_packed
  public :: cholesky_factor_regularized, cholesky_solve

contains

  pure integer function packed_size(n) result(m)
    integer, intent(in) :: n
    m = n * (n + 1) / 2
  end function packed_size

  pure integer function packed_index(i, j) result(k)
    integer, intent(in) :: i, j
    integer :: ii, jj
    ii = min(i, j)
    jj = max(i, j)
    k = ii + jj * (jj - 1) / 2
  end function packed_index

  subroutine sym_packed_matvec(a, x, y, add)
    real(dp), intent(in) :: a(:), x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: add
    integer :: i, j, k, n
    logical :: do_add

    n = size(x)
    do_add = .true.
    if (present(add)) do_add = add
    if (.not. do_add) y = 0.0_dp

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        if (i == j) then
          y(i) = y(i) + a(k) * x(i)
        else
          y(i) = y(i) + a(k) * x(j)
          y(j) = y(j) + a(k) * x(i)
        end if
      end do
    end do
  end subroutine sym_packed_matvec

  subroutine sym_packed_matvec_indexed(a, x, idx, y_local)
    real(dp), intent(in) :: a(:), x(:)
    integer, intent(in) :: idx(:)
    real(dp), intent(out) :: y_local(:)
    integer :: i, j, k, n

    n = size(idx)
    y_local = 0.0_dp
    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        if (i == j) then
          y_local(i) = y_local(i) + a(k) * x(idx(i))
        else
          y_local(i) = y_local(i) + a(k) * x(idx(j))
          y_local(j) = y_local(j) + a(k) * x(idx(i))
        end if
      end do
    end do
  end subroutine sym_packed_matvec_indexed

  subroutine rank_one_update_packed(a, x, scale)
    real(dp), intent(inout) :: a(:)
    real(dp), intent(in) :: x(:), scale
    integer :: i, j, k

    k = 0
    do j = 1, size(x)
      do i = 1, j
        k = k + 1
        a(k) = a(k) + scale * x(i) * x(j)
      end do
    end do
  end subroutine rank_one_update_packed

  subroutine bfgs_inverse_update_packed(h, s, hy, yhy, rho)
    real(dp), intent(inout) :: h(:)
    real(dp), intent(in) :: s(:), hy(:), yhy, rho
    integer :: i, j, k
    real(dp) :: fac

    fac = rho * (rho * yhy + 1.0_dp)
    k = 0
    do j = 1, size(s)
      do i = 1, j
        k = k + 1
        h(k) = h(k) + fac * s(i) * s(j) - rho * (hy(j) * s(i) + hy(i) * s(j))
      end do
    end do
  end subroutine bfgs_inverse_update_packed

  pure subroutine kahan_add(sum_value, compensation, value)
    real(dp), intent(inout) :: sum_value, compensation
    real(dp), intent(in) :: value
    real(dp) :: y, t
    y = value - compensation
    t = sum_value + y
    compensation = (t - sum_value) - y
    sum_value = t
  end subroutine kahan_add

  subroutine packed_to_dense(a, m)
    real(dp), intent(in) :: a(:)
    real(dp), intent(out) :: m(:,:)
    integer :: i, j, k

    m = 0.0_dp
    k = 0
    do j = 1, size(m, 2)
      do i = 1, j
        k = k + 1
        m(i,j) = a(k)
        m(j,i) = a(k)
      end do
    end do
  end subroutine packed_to_dense

  subroutine dense_to_packed(m, a)
    real(dp), intent(in) :: m(:,:)
    real(dp), intent(out) :: a(:)
    integer :: i, j, k
    k = 0
    do j = 1, size(m, 2)
      do i = 1, j
        k = k + 1
        a(k) = m(i,j)
      end do
    end do
  end subroutine dense_to_packed

  subroutine cholesky_factor_regularized(a, l, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: work(:,:)
    real(dp) :: min_diag, shift, eps0
    integer :: i, attempt

    allocate(work(size(a,1), size(a,2)))
    min_diag = minval([(a(i,i), i=1,size(a,1))])
    eps0 = 1.0e-3_dp
    shift = max(0.0_dp, -min_diag)
    ok = .false.

    do attempt = 0, 10
      work = a
      if (attempt > 0 .or. shift > 0.0_dp) then
        do i = 1, size(a,1)
          work(i,i) = work(i,i) + shift + eps0 * 10.0_dp**max(0, attempt-1)
        end do
      end if
      call chol_lower(work, l, ok)
      if (ok) return
    end do

    l = 0.0_dp
    do i = 1, size(a,1)
      l(i,i) = sqrt(abs(a(i,i)) + eps0)
    end do
    ok = .true.
  contains
    subroutine chol_lower(m, c, success)
      real(dp), intent(in) :: m(:,:)
      real(dp), intent(out) :: c(:,:)
      logical, intent(out) :: success
      integer :: ii, jj, kk, n
      real(dp) :: s
      n = size(m,1)
      c = 0.0_dp
      success = .false.
      do ii = 1, n
        do jj = 1, ii
          s = m(ii,jj)
          do kk = 1, jj - 1
            s = s - c(ii,kk) * c(jj,kk)
          end do
          if (ii == jj) then
            if (.not. (s > 0.0_dp)) return
            c(ii,jj) = sqrt(s)
          else
            c(ii,jj) = s / c(jj,jj)
          end if
        end do
      end do
      success = .true.
    end subroutine chol_lower
  end subroutine cholesky_factor_regularized

  subroutine cholesky_solve(l, b, x)
    real(dp), intent(in) :: l(:,:), b(:)
    real(dp), intent(out) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: i, j, n

    n = size(b)
    allocate(y(n))
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
  end subroutine cholesky_solve

end module psqn_linalg
