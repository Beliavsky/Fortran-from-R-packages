! SPDX-License-Identifier: GPL-2.0-or-later
module test_hessian_functions
   use numderiv, only : dp
   implicit none
contains
   function quadratic(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp) :: value
      value = 2.0_dp * z(1) ** 2 + 3.0_dp * z(2) ** 2 + z(3) ** 2 + &
         1.5_dp * z(1) * z(2) - 0.7_dp * z(1) * z(3) + 0.8_dp * z(2) * z(3)
   end function quadratic

   function quadratic_complex(z) result(value)
      complex(dp), intent(in) :: z(:)
      complex(dp) :: value
      value = 2.0_dp * z(1) ** 2 + 3.0_dp * z(2) ** 2 + z(3) ** 2 + &
         1.5_dp * z(1) * z(2) - 0.7_dp * z(1) * z(3) + 0.8_dp * z(2) * z(3)
   end function quadratic_complex

   function vector_model(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp), allocatable :: value(:)
      value = [z(1), z(1) * z(2), z(2) ** 2]
   end function vector_model

end module test_hessian_functions

program test_hessian
   use test_hessian_functions
   use numderiv, only : dp, gend_result, nd_success, hessian, hessian_complex, gend
   implicit none

   real(dp) :: x(3), exact(3, 3)
   real(dp), allocatable :: hess(:, :)
   type(gend_result) :: result
   integer :: status

   x = [0.4_dp, -0.8_dp, 1.2_dp]
   exact = reshape([ &
      4.0_dp, 1.5_dp, -0.7_dp, &
      1.5_dp, 6.0_dp, 0.8_dp, &
      -0.7_dp, 0.8_dp, 2.0_dp], [3, 3])

   call hessian(quadratic, x, hess, status=status)
   call check(status == nd_success, 'Richardson Hessian status')
   call check_close(hess, exact, 2.0e-9_dp, 'Richardson Hessian')

   call hessian_complex(quadratic_complex, x, hess, status=status)
   call check(status == nd_success, 'complex Hessian status')
   call check_close(hess, exact, 2.0e-10_dp, 'complex Hessian')

   call gend(vector_model, [2.0_dp, 3.0_dp], result)
   call check(result%status == nd_success, 'genD status')
   call check(size(result%dmat, 1) == 3 .and. size(result%dmat, 2) == 5, 'genD shape')
   call check_close(result%dmat, reshape([ &
      1.0_dp, 3.0_dp, 0.0_dp, &
      0.0_dp, 2.0_dp, 6.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 2.0_dp], [3, 5]), 2.0e-6_dp, 'genD values')

   print '(a)', 'test_hessian: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:, :), expected(:, :), tolerance
      character(len=*), intent(in) :: label
      call check(maxval(abs(actual - expected)) <= tolerance, label)
   end subroutine check_close

end program test_hessian
