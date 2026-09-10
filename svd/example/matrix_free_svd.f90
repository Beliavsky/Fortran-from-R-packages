! SPDX-License-Identifier: GPL-2.0-or-later
module matrix_free_example_support
   use svd, only : dp
   implicit none
   private
   public :: diagonal_context, mul, tmul

   type :: diagonal_context
      real(dp), allocatable :: diagonal(:)
   end type diagonal_context

contains

   subroutine mul(context, x, y)
      class(*), intent(inout) :: context !! Diagonal-operator state shared by the callbacks.
      real(dp), intent(in) :: x(:) !! Input vector with the diagonal operator's column dimension.
      real(dp), intent(out) :: y(:) !! Forward product vector.

      select type (context)
      type is (diagonal_context)
         y = context%diagonal * x
      class default
         error stop "matrix_free_svd: invalid context"
      end select
   end subroutine mul

   subroutine tmul(context, x, y)
      class(*), intent(inout) :: context !! Diagonal-operator state shared by the callbacks.
      real(dp), intent(in) :: x(:) !! Input vector with the diagonal operator's row dimension.
      real(dp), intent(out) :: y(:) !! Transpose product vector.

      call mul(context, x, y)
   end subroutine tmul

end module matrix_free_example_support

program matrix_free_svd
   use svd
   use matrix_free_example_support
   implicit none

   type(diagonal_context), target :: context
   type(extmat_callbacks) :: callbacks
   type(extmat_operator) :: op
   type(trlan_svd_result) :: fit
   integer :: info

   context%diagonal = [6.0_dp, 4.0_dp, 1.0_dp]
   callbacks%mul => mul
   callbacks%tmul => tmul
   op = make_extmat(callbacks, context, 3, 3, info)
   if (info /= 0) error stop "matrix_free_svd: invalid operator"

   fit = trlan_svd(op, 2)
   if (fit%info /= 0) error stop "matrix_free_svd: solver failed"
   print '(a,*(1x,f8.4))', 'leading singular values:', fit%d
end program matrix_free_svd
