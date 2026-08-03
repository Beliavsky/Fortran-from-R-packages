! SPDX-License-Identifier: GPL-3.0-only
module test_support
  use pa_kinds, only: dp
  implicit none
  private
  public :: assert_true, assert_close, assert_all_close
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual,expected,tolerance
    character(len=*),intent(in)::message
    if(abs(actual-expected)>tolerance) then
      write(*,'(a,2(1x,es16.8))') 'FAIL: '//trim(message),actual,expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_all_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual(:),expected(:),tolerance
    character(len=*),intent(in)::message
    if(size(actual)/=size(expected)) then
      write(*,'(a)') 'FAIL: '//trim(message)//' size mismatch'
      error stop 1
    end if
    if(maxval(abs(actual-expected))>tolerance) then
      write(*,'(a,1x,es16.8)') 'FAIL: '//trim(message),maxval(abs(actual-expected))
      error stop 1
    end if
  end subroutine assert_all_close
end module test_support
