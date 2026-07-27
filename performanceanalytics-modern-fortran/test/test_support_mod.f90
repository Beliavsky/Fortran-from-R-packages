! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module test_support_mod
  use kinds_mod, only: dp
  implicit none
  private
  public :: assert_close, assert_true, assert_vector_close, assert_matrix_close
contains
  subroutine assert_close(actual,expected,tol,label)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::label
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected)))then
      write(*,'(a,2(1x,es24.15))')'FAIL '//trim(label)//':',actual,expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if(.not.condition)then
      write(*,'(a)')'FAIL '//trim(label)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_vector_close(actual,expected,tol,label)
    real(dp),intent(in)::actual(:),expected(:),tol
    character(len=*),intent(in)::label
    call assert_true(size(actual)==size(expected),trim(label)//' size')
    if(size(actual)>0)call assert_close(maxval(abs(actual-expected)),0.0_dp,tol,label)
  end subroutine assert_vector_close
  subroutine assert_matrix_close(actual,expected,tol,label)
    real(dp),intent(in)::actual(:,:),expected(:,:),tol
    character(len=*),intent(in)::label
    call assert_true(all(shape(actual)==shape(expected)),trim(label)//' shape')
    if(size(actual)>0)call assert_close(maxval(abs(actual-expected)),0.0_dp,tol,label)
  end subroutine assert_matrix_close
end module test_support_mod
