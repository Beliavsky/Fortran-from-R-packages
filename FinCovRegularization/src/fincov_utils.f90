! SPDX-License-Identifier: GPL-2.0-only
module fincov_utils
   use fincov_kinds, only : dp
   implicit none
   private
   public :: lowercase, vector_norm2, is_finite_vector, is_finite_matrix
contains
   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

   pure function vector_norm2(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sqrt(sum(x*x))
   end function vector_norm2

   pure logical function is_finite_vector(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:)
      is_finite_vector = all(ieee_is_finite(x))
   end function is_finite_vector

   pure logical function is_finite_matrix(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:,:)
      is_finite_matrix = all(ieee_is_finite(x))
   end function is_finite_matrix
end module fincov_utils
