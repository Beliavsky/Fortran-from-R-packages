module qpxt_test_support
  use quadprog_kinds, only: dp
  implicit none
  private
  public :: assert_true, assert_close, assert_vector_close
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tol, message)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tol) then
      write(*, '(a,es24.16,a,es24.16)') 'FAIL actual=', actual, &
        ' expected=', expected
      write(*, '(a)') trim(message)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual, expected, tol, message)
    real(dp), intent(in) :: actual(:), expected(:), tol
    character(len=*), intent(in) :: message
    call assert_true(size(actual) == size(expected), message // ' size')
    if (size(actual) > 0) then
      if (maxval(abs(actual - expected)) > tol) then
        write(*, '(a,es24.16)') 'FAIL max error=', &
          maxval(abs(actual - expected))
        write(*, '(a)') trim(message)
        error stop 1
      end if
    end if
  end subroutine assert_vector_close
end module qpxt_test_support
