module test_support
  use mvtnorm_kinds, only : dp
  implicit none
contains
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond
    character(len=*),intent(in)::msg
    if(.not.cond) then
      write(*,'(a)') 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol) then
      write(*,'(a,1x,es24.16,1x,es24.16,1x,es12.4)') 'FAIL: '//trim(msg),x,y,abs(x-y)
      error stop 1
    end if
  end subroutine assert_close
end module test_support
