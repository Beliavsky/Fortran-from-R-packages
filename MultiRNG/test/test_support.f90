! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module test_support
  use multirng_kinds, only : dp
  implicit none
  private
  public :: assert_true, assert_close, column_mean, sample_covariance
contains
  subroutine assert_true(cond, msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (.not. cond) then
      write(*, '(a)') 'FAIL: ' // trim(msg)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(x, y, tol, msg)
    real(dp), intent(in) :: x, y, tol
    character(len=*), intent(in) :: msg
    call assert_true(abs(x - y) <= tol, msg)
  end subroutine assert_close

  function column_mean(x) result(m)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable :: m(:)
    integer :: j
    allocate(m(size(x, 2)))
    do j = 1, size(x, 2)
      m(j) = sum(x(:, j)) / real(size(x, 1), dp)
    end do
  end function column_mean

  function sample_covariance(x) result(c)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable :: c(:, :), m(:)
    integer :: i, j, k, n, d
    n = size(x, 1)
    d = size(x, 2)
    allocate(c(d, d))
    m = column_mean(x)
    c = 0.0_dp
    do i = 1, d
      do j = 1, d
        do k = 1, n
          c(i, j) = c(i, j) + (x(k, i) - m(i)) * (x(k, j) - m(j))
        end do
        c(i, j) = c(i, j) / real(n - 1, dp)
      end do
    end do
  end function sample_covariance
end module test_support
