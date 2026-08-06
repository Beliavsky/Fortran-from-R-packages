module test_utils
  use tscopula_kinds, only : dp
  implicit none
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//trim(message);error stop 1;end if
  end subroutine assert_true
  subroutine assert_close(actual,expected,tol,message)
    real(dp),intent(in)::actual,expected,tol;character(len=*),intent(in)::message
    if(abs(actual-expected)>tol)then;write(*,'(a,2es18.8)')'FAIL: '//trim(message)//' actual/expected=',actual,expected;error stop 1;end if
  end subroutine assert_close
  subroutine pass(name)
    character(len=*),intent(in)::name;write(*,'(a)')'PASS: '//trim(name)
  end subroutine pass
end module test_utils
