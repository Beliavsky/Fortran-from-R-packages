module test_support
  use ghyp_kinds, only : dp
  implicit none
  private
  public :: assert_true, assert_close, assert_all_close, finish_test
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
      write(*,'(a,2es18.8)') 'FAIL: '//trim(message)//' actual/expected: ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_all_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    if (size(actual) /= size(expected) .or. maxval(abs(actual - expected)) > tolerance) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_all_close

  subroutine finish_test(name)
    character(len=*), intent(in) :: name
    write(*,'(a)') 'PASS: '//trim(name)
  end subroutine finish_test
end module test_support
