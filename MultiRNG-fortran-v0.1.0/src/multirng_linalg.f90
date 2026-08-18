! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multirng_linalg
  use multirng_kinds, only : dp
  implicit none
  private

  public :: cholesky_lower, invert_spd, is_symmetric_pd

contains

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: l(:, :)
    integer, intent(out) :: info
    integer :: n, i, j, k
    real(dp) :: s

    n = size(a, 1)
    if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
      info = -1
      return
    end if

    l = 0.0_dp
    info = 0
    do i = 1, n
      do j = 1, i
        s = a(i, j)
        do k = 1, j - 1
          s = s - l(i, k) * l(j, k)
        end do
        if (i == j) then
          if (s <= 0.0_dp) then
            info = i
            return
          end if
          l(i, j) = sqrt(s)
        else
          l(i, j) = s / l(j, j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  logical function is_symmetric_pd(a, tol) result(ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: l(:, :)
    real(dp) :: t
    integer :: n, info

    n = size(a, 1)
    if (size(a, 2) /= n) then
      ok = .false.
      return
    end if
    t = 100.0_dp * epsilon(1.0_dp)
    if (present(tol)) t = tol
    if (maxval(abs(a - transpose(a))) > t * max(1.0_dp, maxval(abs(a)))) then
      ok = .false.
      return
    end if
    allocate(l(n, n))
    call cholesky_lower(a, l, info)
    ok = info == 0
  end function is_symmetric_pd

  subroutine invert_spd(a, ainv, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: ainv(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:, :), y(:), x(:)
    integer :: n, i, j, k

    n = size(a, 1)
    if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
      info = -1
      return
    end if

    allocate(l(n, n), y(n), x(n))
    call cholesky_lower(a, l, info)
    if (info /= 0) return

    ainv = 0.0_dp
    do j = 1, n
      y = 0.0_dp
      do i = 1, n
        y(i) = merge(1.0_dp, 0.0_dp, i == j)
        do k = 1, i - 1
          y(i) = y(i) - l(i, k) * y(k)
        end do
        y(i) = y(i) / l(i, i)
      end do

      x = 0.0_dp
      do i = n, 1, -1
        x(i) = y(i)
        do k = i + 1, n
          x(i) = x(i) - l(k, i) * x(k)
        end do
        x(i) = x(i) / l(i, i)
      end do
      ainv(:, j) = x
    end do
    ainv = 0.5_dp * (ainv + transpose(ainv))
  end subroutine invert_spd

end module multirng_linalg
