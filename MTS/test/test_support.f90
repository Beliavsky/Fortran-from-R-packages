! SPDX-License-Identifier: Artistic-2.0
module test_support
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mts_kinds, only : dp
   implicit none
   private
   public :: assert_true, assert_close, assert_matrix_close, assert_finite
contains
   subroutine assert_true(condition,message)
      logical,intent(in)::condition
      character(len=*),intent(in)::message
      if(.not.condition)then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual,expected,tolerance,message)
      real(dp),intent(in)::actual,expected,tolerance
      character(len=*),intent(in)::message
      if(.not.ieee_is_finite(actual).or.abs(actual-expected)>tolerance)then
         write(*,'(a,2(1x,es16.8))') 'FAIL: '//trim(message)//' actual expected:',actual,expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_matrix_close(actual,expected,tolerance,message)
      real(dp),intent(in)::actual(:,:),expected(:,:),tolerance
      character(len=*),intent(in)::message
      if(any(shape(actual)/=shape(expected)).or.any(.not.ieee_is_finite(actual)).or.maxval(abs(actual-expected))>tolerance)then
         write(*,'(a,1x,es16.8)') 'FAIL: '//trim(message)//' max error:',maxval(abs(actual-expected))
         error stop 1
      end if
   end subroutine assert_matrix_close

   subroutine assert_finite(x,message)
      real(dp),intent(in)::x(:)
      character(len=*),intent(in)::message
      if(any(.not.ieee_is_finite(x)))then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine assert_finite
end module test_support
