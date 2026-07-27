! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module test_support
  use tsdyn_kinds, only: dp
  implicit none
  private
  public :: assert_true, assert_close, assert_finite, assert_all_finite
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(x,y,tol,message)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::message
    call assert_true(abs(x-y)<=tol,trim(message))
  end subroutine assert_close
  subroutine assert_finite(x,message)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp),intent(in)::x
    character(len=*),intent(in)::message
    call assert_true(ieee_is_finite(x),message)
  end subroutine assert_finite
  subroutine assert_all_finite(x,message)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp),intent(in)::x(..)
    character(len=*),intent(in)::message
    select rank(x)
    rank(1);call assert_true(all(ieee_is_finite(x)),message)
    rank(2);call assert_true(all(ieee_is_finite(x)),message)
    rank(3);call assert_true(all(ieee_is_finite(x)),message)
    rank default;call assert_true(.false.,message)
    end select
  end subroutine assert_all_finite
end module test_support
