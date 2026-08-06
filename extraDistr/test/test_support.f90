module test_support
  use extra_distr_kinds, only : dp
  implicit none
  private
  public :: assert_close, assert_true, assert_int, finish_tests
  integer :: failures = 0
contains
  subroutine assert_close(actual,expected,tol,message)
    real(dp),intent(in)::actual,expected,tol
    character(*),intent(in)::message
    if(.not.(abs(actual-expected)<=tol))then
      failures=failures+1
      write(*,'(a,2(1x,es24.15))') 'FAIL '//trim(message)//':',actual,expected
    end if
  end subroutine assert_close
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(*),intent(in)::message
    if(.not.condition)then
      failures=failures+1
      write(*,'(a)') 'FAIL '//trim(message)
    end if
  end subroutine assert_true
  subroutine assert_int(actual,expected,message)
    integer,intent(in)::actual,expected
    character(*),intent(in)::message
    if(actual/=expected)then
      failures=failures+1
      write(*,'(a,2(1x,i0))') 'FAIL '//trim(message)//':',actual,expected
    end if
  end subroutine assert_int
  subroutine finish_tests(name)
    character(*),intent(in)::name
    if(failures>0)then
      write(*,'(a,i0)') trim(name)//' failures: ',failures
      error stop 1
    end if
    write(*,'(a)') trim(name)//': PASS'
  end subroutine finish_tests
end module test_support
