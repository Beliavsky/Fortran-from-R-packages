! SPDX-License-Identifier: GPL-3.0-only
program test_sparse
   use matrix, only : dp, csr_matrix, csc_matrix, csr_from_triplet, csr_from_dense, csr_to_dense, &
      csc_from_csr, csr_from_csc, csr_transpose, csr_matvec, csr_add, csr_multiply, &
      csr_kronecker, csr_permute, csr_is_symmetric, csr_row_sums, csr_col_sums, &
      csr_crossprod, csr_lu_factorize, csr_lu_solve, sparse_lu_factor, eye, matrix_success
   implicit none
   type(csr_matrix) :: a, b, c, at, back, cp
   type(csc_matrix) :: ac
   type(sparse_lu_factor) :: fac
   real(dp), allocatable :: dense(:,:), y(:), rhs(:,:), x(:,:)
   integer :: info

   call csr_from_triplet(3, 3, [1, 1, 2, 3, 3], [1, 1, 2, 1, 3], &
      [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], a, info)
   call check(info == matrix_success .and. a%nnz() == 4 .and. a%valid(), 'triplet')
   dense = csr_to_dense(a)
   call check_close(dense, reshape([3.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 6.0_dp], [3, 3]), 1.0e-14_dp, 'to dense')

   call csr_matvec(a, [1.0_dp, 2.0_dp, 3.0_dp], y, info)
   call check_close_vector(y, [3.0_dp, 8.0_dp, 23.0_dp], 1.0e-14_dp, 'matvec')
   call csr_transpose(a, at)
   call csc_from_csr(a, ac)
   call csr_from_csc(ac, back)
   call check_close(csr_to_dense(at), transpose(dense), 1.0e-14_dp, 'transpose')
   call check_close(csr_to_dense(back), dense, 1.0e-14_dp, 'csc roundtrip')

   call csr_add(a, at, b, info)
   call check_close(csr_to_dense(b), dense + transpose(dense), 1.0e-14_dp, 'add')
   call csr_multiply(a, at, c, info)
   call check_close(csr_to_dense(c), matmul(dense, transpose(dense)), 1.0e-13_dp, 'multiply')
   call csr_crossprod(a, cp, info)
   call check_close(csr_to_dense(cp), matmul(transpose(dense), dense), 1.0e-13_dp, 'crossprod')

   call csr_from_dense(eye(2), b)
   call csr_kronecker(b, a, c, info)
   call check(size(csr_to_dense(c), 1) == 6 .and. c%nnz() == 8, 'sparse kronecker')
   call csr_permute(a, b, info, row_perm=[3, 2, 1], col_perm=[3, 2, 1])
   call check_close(csr_to_dense(b), dense([3, 2, 1], [3, 2, 1]), 1.0e-14_dp, 'sparse permute')
   call check(.not. csr_is_symmetric(a), 'nonsymmetric')
   call check_close_vector(csr_row_sums(a), [3.0_dp, 4.0_dp, 11.0_dp], 1.0e-14_dp, 'row sums')
   call check_close_vector(csr_col_sums(a), [8.0_dp, 4.0_dp, 6.0_dp], 1.0e-14_dp, 'col sums')

   call csr_from_dense(reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2]), b)
   call csr_lu_factorize(b, fac, info)
   rhs = reshape([1.0_dp, 2.0_dp], [2, 1])
   call csr_lu_solve(fac, rhs, x, info)
   call check_close(matmul(csr_to_dense(b), x), rhs, 1.0e-12_dp, 'sparse solve')

   print '(a)', 'test_sparse: PASS'
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

   subroutine check_close_vector(actual, expected, tol, name)
      real(dp), intent(in) :: actual(:), expected(:)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: name
      call check(size(actual) == size(expected) .and. maxval(abs(actual - expected)) <= tol, name)
   end subroutine check_close_vector
end program test_sparse
