! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module test_support
  use rq_kinds, only: dp
  implicit none
  private
  public :: assert_close, assert_true
contains
  subroutine assert_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual,expected,tolerance
    character(len=*),intent(in)::message
    if(abs(actual-expected)>tolerance) then
      write(*,'(a)') 'FAIL: '//trim(message)
      write(*,'(a,es24.16)') ' actual:   ',actual
      write(*,'(a,es24.16)') ' expected: ',expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true
end module test_support
