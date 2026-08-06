! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_test_support
  use waveslim_kinds, only : dp
  implicit none
  private
  public :: assert_true, assert_close_scalar, assert_close_array
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close_scalar(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual-expected) > tolerance) then
      write(*,'(a,3(1x,es14.6))') 'FAIL: '//trim(message), actual, expected, tolerance
      error stop 1
    end if
  end subroutine assert_close_scalar

  subroutine assert_close_array(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    real(dp) :: error_value
    if (size(actual) /= size(expected)) then
      write(*,'(a)') 'FAIL: '//trim(message)//' shape mismatch'
      error stop 1
    end if
    error_value = maxval(abs(actual-expected))
    if (error_value > tolerance) then
      write(*,'(a,2(1x,es14.6))') 'FAIL: '//trim(message), error_value, tolerance
      error stop 1
    end if
  end subroutine assert_close_array
end module waveslim_test_support
