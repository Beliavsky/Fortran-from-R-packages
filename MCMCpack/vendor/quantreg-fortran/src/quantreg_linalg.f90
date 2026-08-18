! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_linalg
  use quantreg_kinds, only : dp
  implicit none
  private
  public :: spd_solve, symmetric_inverse, solve_linear
contains

  subroutine spd_solve(a, b, info)
    real(dp), intent(inout) :: a(:,:)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info
    integer :: i, j, k, n
    real(dp) :: s

    n = size(b)
    info = 0
    if (size(a,1) /= n .or. size(a,2) /= n) then
      info = -1
      return
    end if

    do j = 1, n
      s = a(j,j)
      do k = 1, j - 1
        s = s - a(j,k) * a(j,k)
      end do
      if (s <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        info = j
        return
      end if
      a(j,j) = sqrt(s)
      do i = j + 1, n
        s = a(i,j)
        do k = 1, j - 1
          s = s - a(i,k) * a(j,k)
        end do
        a(i,j) = s / a(j,j)
      end do
    end do

    do i = 1, n
      s = b(i)
      do k = 1, i - 1
        s = s - a(i,k) * b(k)
      end do
      b(i) = s / a(i,i)
    end do
    do i = n, 1, -1
      s = b(i)
      do k = i + 1, n
        s = s - a(k,i) * b(k)
      end do
      b(i) = s / a(i,i)
    end do
  end subroutine spd_solve

  subroutine solve_linear(a, b, info)
    real(dp), intent(inout) :: a(:,:)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info
    integer :: i, j, k, n, piv
    real(dp) :: best, fac, tmp

    n = size(b)
    info = 0
    do k = 1, n - 1
      piv = k
      best = abs(a(k,k))
      do i = k + 1, n
        if (abs(a(i,k)) > best) then
          best = abs(a(i,k))
          piv = i
        end if
      end do
      if (best <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        info = k
        return
      end if
      if (piv /= k) then
        do j = k, n
          tmp = a(k,j); a(k,j) = a(piv,j); a(piv,j) = tmp
        end do
        tmp = b(k); b(k) = b(piv); b(piv) = tmp
      end if
      do i = k + 1, n
        fac = a(i,k) / a(k,k)
        a(i,k) = 0.0_dp
        a(i,k+1:n) = a(i,k+1:n) - fac * a(k,k+1:n)
        b(i) = b(i) - fac * b(k)
      end do
    end do
    if (abs(a(n,n)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
      info = n
      return
    end if
    do i = n, 1, -1
      tmp = b(i)
      if (i < n) tmp = tmp - dot_product(a(i,i+1:n), b(i+1:n))
      b(i) = tmp / a(i,i)
    end do
  end subroutine solve_linear

  subroutine symmetric_inverse(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: work(:,:), rhs(:)
    integer :: j, n

    n = size(a,1)
    info = 0
    do j = 1, n
      allocate(work(n,n), rhs(n))
      work = a
      rhs = 0.0_dp
      rhs(j) = 1.0_dp
      call spd_solve(work, rhs, info)
      if (info /= 0) return
      ainv(:,j) = rhs
      deallocate(work, rhs)
    end do
    ainv = 0.5_dp * (ainv + transpose(ainv))
  end subroutine symmetric_inverse

end module quantreg_linalg
