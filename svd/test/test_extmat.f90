! SPDX-License-Identifier: GPL-2.0-or-later
module test_extmat_support
   use svd, only : dp
   implicit none
   private
   public :: matrix_context, forward_product, transpose_product

   type :: matrix_context
      real(dp), allocatable :: a(:, :)
   end type matrix_context

contains

   subroutine forward_product(context, x, y)
      class(*), intent(inout) :: context !! Test context storing the dense reference matrix.
      real(dp), intent(in) :: x(:) !! Input vector with length equal to the reference matrix column count.
      real(dp), intent(out) :: y(:) !! Product vector with length equal to the reference matrix row count.

      select type (context)
      type is (matrix_context)
         y = matmul(context%a, x)
      class default
         error stop "test_extmat: unexpected context type"
      end select
   end subroutine forward_product

   subroutine transpose_product(context, x, y)
      class(*), intent(inout) :: context !! Test context storing the dense reference matrix.
      real(dp), intent(in) :: x(:) !! Input vector with length equal to the reference matrix row count.
      real(dp), intent(out) :: y(:) !! Transpose-product vector with length equal to the reference matrix column count.

      select type (context)
      type is (matrix_context)
         y = matmul(transpose(context%a), x)
      class default
         error stop "test_extmat: unexpected context type"
      end select
   end subroutine transpose_product

end module test_extmat_support

program test_extmat
   use rspectra, only : dense_operator, make_dense_operator
   use svd
   use test_extmat_support
   implicit none

   type(matrix_context), target :: context
   type(extmat_callbacks) :: callbacks
   type(extmat_operator) :: op
   type(dense_operator) :: dense
   real(dp), allocatable :: y(:), materialized(:, :), values(:)
   real(dp) :: v2(2), v3(3)
   integer :: info

   context%a = reshape([1.0_dp, 3.0_dp, 5.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], [3, 2])
   callbacks%mul => forward_product
   callbacks%tmul => transpose_product
   op = make_extmat(callbacks, context, 3, 2, info)

   call check(info == 0, "constructor status")
   call check(is_extmat(op), "is_extmat true")
   call check(extmat_nrow(op) == 3, "row count")
   call check(extmat_ncol(op) == 2, "column count")

   dense = make_dense_operator(context%a)
   call check(.not. is_extmat(dense), "dense operator is not extmat")

   v2 = [2.0_dp, -1.0_dp]
   call ematmul(op, v2, y, info=info)
   call check(info == 0, "forward status")
   call check(maxval(abs(y - matmul(context%a, v2))) < 1.0e-13_dp, "forward product")

   v3 = [1.0_dp, -2.0_dp, 0.5_dp]
   call ematmul(op, v3, y, transposed=.true., info=info)
   call check(info == 0, "transpose status")
   call check(maxval(abs(y - matmul(transpose(context%a), v3))) < 1.0e-13_dp, "transpose product")

   call materialize_extmat(op, materialized, info)
   call check(info == 0, "materialization status")
   call check(maxval(abs(materialized - context%a)) < 1.0e-13_dp, "materialization")

   call extmat_vector(op, values, info)
   call check(info == 0, "vectorization status")
   call check(maxval(abs(values - reshape(context%a, [size(context%a)]))) < 1.0e-13_dp, "vectorization")

   print '(a)', 'test_extmat: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must be true.
      character(len=*), intent(in) :: label !! Short assertion label printed on failure.

      if (.not. condition) then
         print '(a,1x,a)', 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_extmat
