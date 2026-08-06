module test_support
   use fgarch_kinds, only : dp
   implicit none
   private
   public :: assert_true, assert_close, assert_all_close
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected) > tolerance*(1.0_dp+abs(expected))) then
         write(*,'(a,2es24.14)') 'FAIL: '//trim(message)//' actual/expected: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_all_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: message
      integer :: i
      call assert_true(size(actual) == size(expected),trim(message)//' size')
      do i = 1, size(actual)
         call assert_close(actual(i),expected(i),tolerance,trim(message))
      end do
   end subroutine assert_all_close
end module test_support
