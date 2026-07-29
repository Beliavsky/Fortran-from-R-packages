! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

module test_support
   use fracdiff_kinds, only : dp
   implicit none
   private
   public :: assert_true, assert_close, assert_vector_close, assert_matrix_symmetric

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') "FAIL: "//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected) > tolerance*max(1.0_dp,abs(expected))) then
         write(*,'(a)') "FAIL: "//trim(message)
         write(*,'(a,es24.15)') " actual:   ", actual
         write(*,'(a,es24.15)') " expected: ", expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: message
      integer :: i
      call assert_true(size(actual)==size(expected), trim(message)//" size")
      do i=1,size(actual)
         if (abs(actual(i)-expected(i)) > tolerance*max(1.0_dp,abs(expected(i)))) then
            write(*,'(a,i0)') "FAIL: "//trim(message)//" at index ",i
            write(*,'(a,es24.15)') " actual:   ", actual(i)
            write(*,'(a,es24.15)') " expected: ", expected(i)
            error stop 1
         end if
      end do
   end subroutine assert_vector_close

   subroutine assert_matrix_symmetric(matrix,tolerance,message)
      real(dp),intent(in)::matrix(:,:),tolerance
      character(len=*),intent(in)::message
      call assert_true(size(matrix,1)==size(matrix,2),trim(message)//" square")
      call assert_true(maxval(abs(matrix-transpose(matrix))) <= tolerance, message)
   end subroutine assert_matrix_symmetric
end module test_support
