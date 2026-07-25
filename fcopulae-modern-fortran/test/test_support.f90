! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module test_support
  use fcopulae_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: assert_true, assert_close, assert_finite, assert_all_finite
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)')'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,tolerance,message)
    real(dp),intent(in)::actual,expected,tolerance
    character(len=*),intent(in)::message
    if(.not.ieee_is_finite(actual).or.abs(actual-expected)>tolerance)then
      write(*,'(a,2(es24.15,1x),a,es12.4)')'FAIL: '//trim(message)//' actual/expected ',actual,expected,' tolerance ',tolerance
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_finite(value,message)
    real(dp),intent(in)::value
    character(len=*),intent(in)::message
    call assert_true(ieee_is_finite(value),message)
  end subroutine assert_finite

  subroutine assert_all_finite(values,message)
    real(dp),intent(in)::values(:,:)
    character(len=*),intent(in)::message
    call assert_true(all(ieee_is_finite(values)),message)
  end subroutine assert_all_finite
end module test_support
