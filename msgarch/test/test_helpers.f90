! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module test_helpers
   use msgarch_kinds, only : dp
   implicit none
   private
   public :: assert_true, assert_close, assert_all_finite
contains
   subroutine assert_true(condition,message)
      logical,intent(in)::condition
      character(len=*),intent(in)::message
      if(.not.condition)then
         write(*,'(a)')'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual,expected,tolerance,message)
      real(dp),intent(in)::actual,expected,tolerance
      character(len=*),intent(in)::message
      if(abs(actual-expected)>tolerance)then
         write(*,'(a,2(1x,es16.8))')'FAILED: '//trim(message),actual,expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_all_finite(x,message)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp),intent(in)::x(..)
      character(len=*),intent(in)::message
      logical::ok
      ok=.false.
      select rank(x)
      rank(1);ok=all(ieee_is_finite(x))
      rank(2);ok=all(ieee_is_finite(x))
      rank(3);ok=all(ieee_is_finite(x))
      end select
      call assert_true(ok,message)
   end subroutine assert_all_finite
end module test_helpers
