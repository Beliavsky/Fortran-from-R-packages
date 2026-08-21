! SPDX-License-Identifier: GPL-3.0-only
program test_functions
   use matrix, only : dp, matrix_exponential, matrix_power, matrix_sqrt_sym, &
      near_positive_definite, condition_number_1, reciprocal_condition_1, &
      symmetric_eigen, eye, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:), values(:), vectors(:,:)
   real(dp) :: cond, rcond
   integer :: info

   allocate(a(2, 2), source=0.0_dp)
   a(1, 1) = log(2.0_dp)
   a(2, 2) = log(3.0_dp)
   call matrix_exponential(a, b, info)
   call check(info == matrix_success, 'exponential status')
   call check_close(b, reshape([2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2]), 1.0e-12_dp, 'exponential')

   a = reshape([2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2])
   call matrix_power(a, 3, b, info)
   call check_close(b, reshape([8.0_dp, 0.0_dp, 0.0_dp, 27.0_dp], [2, 2]), 1.0e-12_dp, 'power')
   call matrix_power(a, -1, b, info)
   call check_close(matmul(a, b), eye(2), 1.0e-12_dp, 'negative power')

   call matrix_sqrt_sym(a, b, info)
   call check(info == matrix_success, 'sqrt status')
   call check_close(matmul(b, b), a, 1.0e-11_dp, 'sqrt')

   a = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
   call near_positive_definite(a, b, info, corr=.true., tol=1.0e-10_dp)
   call check(info == matrix_success, 'near pd status')
   call symmetric_eigen(b, values, vectors, info)
   call check(minval(values) > 0.0_dp, 'near pd eigenvalues')
   call check(abs(b(1, 1) - 1.0_dp) < 1.0e-12_dp .and. abs(b(2, 2) - 1.0_dp) < 1.0e-12_dp, &
      'near pd correlation diagonal')

   a = reshape([1.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2])
   call condition_number_1(a, cond, info)
   call reciprocal_condition_1(a, rcond, info)
   call check(abs(cond - 4.0_dp) < 1.0e-12_dp, 'condition number')
   call check(abs(rcond - 0.25_dp) < 1.0e-12_dp, 'reciprocal condition')

   print '(a)', 'test_functions: PASS'
contains
   subroutine check(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not. condition) then
         print '(a)', 'FAIL: ' // name
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tol, name)
      real(dp), intent(in) :: actual(:,:), expected(:,:)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: name
      call check(all(shape(actual) == shape(expected)) .and. maxval(abs(actual - expected)) <= tol, name)
   end subroutine check_close
end program test_functions
