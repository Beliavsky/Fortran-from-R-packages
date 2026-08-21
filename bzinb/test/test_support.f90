module test_support
  use bzinb_kinds, only : dp
  implicit none
contains
  subroutine assert_close(name,a,b,tol,fail)
    character(len=*),intent(in)::name
    real(dp),intent(in)::a,b,tol
    integer,intent(inout)::fail
    if(abs(a-b)>tol)then
      print '(a,2(1x,es14.6),a,1x,es14.6)',trim(name)//' FAIL:',a,b,'tol',tol
      fail=fail+1
    end if
  end subroutine
  subroutine assert_true(name,cond,fail)
    character(len=*),intent(in)::name
    logical,intent(in)::cond
    integer,intent(inout)::fail
    if(.not.cond)then;print '(a)',trim(name)//' FAIL';fail=fail+1;end if
  end subroutine
end module test_support
