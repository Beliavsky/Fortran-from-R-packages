! SPDX-License-Identifier: GPL-3.0-only
module test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mass, only : dp
  implicit none
  private
  public :: assert_true, assert_close, assert_all_finite
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      write(*,'(a,2(1x,es16.8))') 'FAIL: '//trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_all_finite(values, message)
    real(dp), intent(in) :: values(:)
    character(len=*), intent(in) :: message
    call assert_true(all(ieee_is_finite(values)), message)
  end subroutine assert_all_finite
end module test_support
