! SPDX-License-Identifier: GPL-2.0-only
module test_support_block_top
   use streg, only : dp
   implicit none
   private
   public :: assert_true, assert_close, assert_all_finite
contains
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if(.not.condition)then
         write(*,'(a)')'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true
   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual,expected,tolerance
      character(len=*), intent(in) :: message
      call assert_true(abs(actual-expected)<=tolerance*(1.0_dp+abs(expected)),message)
   end subroutine assert_close
   subroutine assert_all_finite(x,message)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:)
      character(len=*), intent(in) :: message
      call assert_true(all(ieee_is_finite(x)),message)
   end subroutine assert_all_finite
end module test_support_block_top
