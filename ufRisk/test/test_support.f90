! SPDX-License-Identifier: GPL-3.0-only
module test_support
   use kind_mod, only : dp
   implicit none
   private
   public :: assert_true, assert_close, assert_vector_finite, finish_tests
   integer :: failures = 0
contains
   subroutine assert_true(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) then
         failures = failures+1
         write(*,'(a)') 'FAIL: '//trim(label)
      end if
   end subroutine assert_true

   subroutine assert_close(actual,expected,tolerance,label)
      real(dp), intent(in) :: actual,expected,tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance*max(1.0_dp,abs(expected))) then
         failures = failures+1
         write(*,'(a,2(1x,es24.15))') 'FAIL: '//trim(label),actual,expected
      end if
   end subroutine assert_close

   subroutine assert_vector_finite(x,label)
      real(dp), intent(in) :: x(:)
      character(len=*), intent(in) :: label
      call assert_true(all(abs(x)<=huge(1.0_dp)),label)
   end subroutine assert_vector_finite

   subroutine finish_tests(name)
      character(len=*), intent(in) :: name
      if (failures > 0) then
         write(*,'(a,i0)') trim(name)//' failures: ',failures
         error stop 1
      end if
      write(*,'(a)') trim(name)//': all tests passed'
   end subroutine finish_tests
end module test_support
