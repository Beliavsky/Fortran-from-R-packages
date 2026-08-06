! SPDX-License-Identifier: GPL-3.0-only
program test_constructors
   use matrix, only : dp, toeplitz_matrix, hilbert_matrix, permutation_matrix, &
      sparse_identity, sparse_diagonal, random_sparse_matrix, companion_matrix, &
      csr_matrix, csr_to_dense, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:)
   type(csr_matrix) :: s
   integer :: info

   a = toeplitz_matrix([1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 4.0_dp, 5.0_dp], info)
   call check(info == matrix_success, 'toeplitz status')
   call check_close(a, reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, &
      5.0_dp, 4.0_dp, 1.0_dp], [3, 3]), 1.0e-14_dp, 'toeplitz')
   a = hilbert_matrix(2)
   call check_close(a, reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp / 3.0_dp], [2, 2]), 1.0e-14_dp, 'hilbert')
   a = permutation_matrix([2, 3, 1], info)
   call check(info == matrix_success .and. abs(sum(a) - 3.0_dp) < 1.0e-14_dp, 'permutation')

   call sparse_identity(3, s)
   call check_close(csr_to_dense(s), reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3]), 1.0e-14_dp, 'sparse identity')
   call sparse_diagonal([2.0_dp, 3.0_dp], s, nrow=3, ncol=2)
   call check_close(csr_to_dense(s), reshape([2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp], &
      [3, 2]), 1.0e-14_dp, 'sparse diagonal')

   call random_sparse_matrix(5, 5, 0.3_dp, s, info, symmetric=.true.)
   b = csr_to_dense(s)
   call check(info == matrix_success .and. maxval(abs(b - transpose(b))) < 1.0e-14_dp, 'random symmetric')

   a = companion_matrix([1.0_dp, -2.0_dp, 3.0_dp])
   call check_close(a, reshape([-1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, &
      -3.0_dp, 0.0_dp, 0.0_dp], [3, 3]), 1.0e-14_dp, 'companion')

   print '(a)', 'test_constructors: PASS'
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
end program test_constructors
