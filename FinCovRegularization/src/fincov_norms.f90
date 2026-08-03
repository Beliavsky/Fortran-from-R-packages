! SPDX-License-Identifier: GPL-2.0-only
module fincov_norms
   use fincov_kinds, only : dp
   use fincov_linalg, only : frobenius_norm_squared, spectral_norm_squared
   implicit none
   private

   public :: f_norm2, o_norm2
contains
   pure function f_norm2(matrix) result(value)
      real(dp), intent(in) :: matrix(:,:)
      real(dp) :: value
      value = frobenius_norm_squared(matrix)
   end function f_norm2

   function o_norm2(matrix, status) result(value)
      real(dp), intent(in) :: matrix(:,:)
      integer, intent(out), optional :: status
      real(dp) :: value
      value = spectral_norm_squared(matrix, status)
   end function o_norm2
end module fincov_norms
