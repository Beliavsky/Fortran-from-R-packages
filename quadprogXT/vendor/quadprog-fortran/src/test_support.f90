! SPDX-License-Identifier: GPL-2.0-or-later
module quadprog_test_support
  use quadprog_kinds, only: dp
  implicit none
  private
  public :: check, check_close, check_vector_close
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    call check(abs(actual - expected) <= tolerance, message)
  end subroutine check_close

  subroutine check_vector_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    call check(size(actual) == size(expected), message // ' size')
    call check(maxval(abs(actual - expected)) <= tolerance, message)
  end subroutine check_vector_close
end module quadprog_test_support
