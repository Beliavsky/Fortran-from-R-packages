! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module test_support
  use fbasics_kinds, only: dp
  implicit none
  private
  public :: assert_close, assert_true, assert_all_finite
contains
  subroutine assert_close(actual,expected,tol,message)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::message
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected)))then
      write(*,'(a,2(1x,es24.15))')'FAIL '//trim(message)//':',actual,expected
      error stop 1
    end if
  end subroutine
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL '//trim(message);error stop 1;end if
  end subroutine
  subroutine assert_all_finite(x,message)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp),intent(in)::x(:)
    character(len=*),intent(in)::message
    if(.not.all(ieee_is_finite(x)))then;write(*,'(a)')'FAIL '//trim(message);error stop 1;end if
  end subroutine
end module test_support
