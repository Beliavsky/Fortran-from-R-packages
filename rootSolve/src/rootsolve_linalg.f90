! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_linalg
  use rootsolve_kinds, only : dp
  implicit none
  private
  public :: solve_linear, norm_inf, vector_norm2
contains

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: piv, fac, tmp
    integer :: n, i, j, k, p

    n = size(b)
    info = 0
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
      info = -1
      x = 0.0_dp
      return
    end if
    allocate(aa(n,n), bb(n))
    aa = a
    bb = b

    do k = 1, n-1
      p = k
      piv = abs(aa(k,k))
      do i = k+1, n
        if (abs(aa(i,k)) > piv) then
          p = i
          piv = abs(aa(i,k))
        end if
      end do
      if (piv <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa)))) then
        info = k
        x = 0.0_dp
        return
      end if
      if (p /= k) then
        do j = k, n
          tmp = aa(k,j)
          aa(k,j) = aa(p,j)
          aa(p,j) = tmp
        end do
        tmp = bb(k)
        bb(k) = bb(p)
        bb(p) = tmp
      end if
      do i = k+1, n
        fac = aa(i,k) / aa(k,k)
        aa(i,k) = 0.0_dp
        do j = k+1, n
          aa(i,j) = aa(i,j) - fac * aa(k,j)
        end do
        bb(i) = bb(i) - fac * bb(k)
      end do
    end do

    if (abs(aa(n,n)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa)))) then
      info = n
      x = 0.0_dp
      return
    end if

    x(n) = bb(n) / aa(n,n)
    do i = n-1, 1, -1
      x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
    end do
  end subroutine solve_linear

  pure real(dp) function norm_inf(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = 0.0_dp
    else
      v = maxval(abs(x))
    end if
  end function norm_inf

  pure real(dp) function vector_norm2(x) result(v)
    real(dp), intent(in) :: x(:)
    v = sqrt(dot_product(x, x))
  end function vector_norm2
end module rootsolve_linalg
