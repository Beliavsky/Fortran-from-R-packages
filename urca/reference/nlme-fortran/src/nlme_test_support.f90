! SPDX-License-Identifier: GPL-2.0-or-later
module test_support
  use nlme_kinds, only : dp
  implicit none
  private
  public :: check, check_close, check_vector_close
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)')'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine check
  subroutine check_close(actual,expected,tol,message)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::message
    call check(abs(actual-expected)<=tol*max(1.0_dp,abs(expected)),message)
  end subroutine check_close
  subroutine check_vector_close(actual,expected,tol,message)
    real(dp),intent(in)::actual(:),expected(:),tol
    character(len=*),intent(in)::message
    call check(size(actual)==size(expected),trim(message)//' size')
    call check(maxval(abs(actual-expected))<=tol*max(1.0_dp,maxval(abs(expected))),message)
  end subroutine check_vector_close
end module test_support
