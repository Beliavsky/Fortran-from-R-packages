! SPDX-License-Identifier: GPL-3.0-only
program test_dense
   use matrix, only : dp, eye, diag_matrix, band_matrix, kronecker_product, khatri_rao, &
      symmpart, skewpart, row_sums, col_sums, pack_triangular, unpack_triangular, &
      dense_permute, trace_matrix, frobenius_norm, crossprod, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), x(:), rs(:), cs(:)
   integer :: info

   a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [2, 2])
   call check_close(eye(3), reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3]), 1.0e-14_dp, 'eye')
   b = diag_matrix([2.0_dp, 3.0_dp])
   call check_close(b, reshape([2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2]), 1.0e-14_dp, 'diag')
   b = band_matrix(a, 0, 1)
   call check_close(b, reshape([1.0_dp, 0.0_dp, 2.0_dp, 4.0_dp], [2, 2]), 1.0e-14_dp, 'band')
   c = kronecker_product(a, reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))
   call check(size(c, 1) == 4 .and. size(c, 2) == 4, 'kronecker shape')
   call check_close(c(1:2, 3:4), 2.0_dp * eye(2), 1.0e-14_dp, 'kronecker block')

   c = khatri_rao(reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2]), &
      reshape([5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp], [2, 2]), info)
   call check(info == matrix_success, 'khatri rao status')
   call check_close(c(:, 1), [5.0_dp, 6.0_dp, 10.0_dp, 12.0_dp], 1.0e-14_dp, 'khatri rao')

   b = symmpart(a, info)
   call check(info == matrix_success, 'symmpart status')
   d = skewpart(a, info)
   call check_close(b + d, a, 1.0e-14_dp, 'sym plus skew')
   rs = row_sums(a)
   cs = col_sums(a)
   call check_close(rs, [3.0_dp, 7.0_dp], 1.0e-14_dp, 'row sums')
   call check_close(cs, [4.0_dp, 6.0_dp], 1.0e-14_dp, 'col sums')

   x = pack_triangular(a, upper=.true., info=info)
   b = unpack_triangular(x, 2, upper=.true., info=info)
   call check_close(b, reshape([1.0_dp, 0.0_dp, 2.0_dp, 4.0_dp], [2, 2]), 1.0e-14_dp, 'pack')
   b = dense_permute(a, row_perm=[2, 1], col_perm=[2, 1], info=info)
   call check_close(b, reshape([4.0_dp, 2.0_dp, 3.0_dp, 1.0_dp], [2, 2]), 1.0e-14_dp, 'permute')
   call check(abs(trace_matrix(a) - 5.0_dp) < 1.0e-14_dp, 'trace')
   call check(abs(frobenius_norm(a) - sqrt(30.0_dp)) < 1.0e-14_dp, 'frobenius')
   c = crossprod(a)
   call check_close(c, matmul(transpose(a), a), 1.0e-14_dp, 'crossprod')

   print '(a)', 'test_dense: PASS'
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
      real(dp), intent(in) :: actual(..), expected(..)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: name
      select rank(actual)
      rank(1)
         select rank(expected)
         rank(1)
            call check(size(actual) == size(expected) .and. &
               maxval(abs(actual - expected)) <= tol, name)
         rank default
            call check(.false., name)
         end select
      rank(2)
         select rank(expected)
         rank(2)
            call check(all(shape(actual) == shape(expected)) .and. &
               maxval(abs(actual - expected)) <= tol, name)
         rank default
            call check(.false., name)
         end select
      rank default
         call check(.false., name)
      end select
   end subroutine check_close
end program test_dense
