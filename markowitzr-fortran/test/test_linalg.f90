! SPDX-License-Identifier: LGPL-3.0-or-later
program test_linalg
   use markowitzr, only: dp, symmetric_vech, symmetric_ivech
   use markowitzr, only: duplication_matrix, kronecker_product
   implicit none
   real(dp) :: a(3,3), vectorized(9), expected(9)
   real(dp), allocatable :: packed(:), restored(:, :), d(:, :), kron(:, :)
   integer :: status

   a = reshape([2.0_dp,1.0_dp,3.0_dp, &
                1.0_dp,4.0_dp,5.0_dp, &
                3.0_dp,5.0_dp,6.0_dp],[3,3])
   packed = symmetric_vech(a)
   call assert_close(packed,[2.0_dp,1.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],1.0e-14_dp)
   restored = symmetric_ivech(packed,status)
   if (status /= 0) error stop 1
   call assert_close(reshape(restored,[9]),reshape(a,[9]),1.0e-14_dp)

   d = duplication_matrix(3)
   vectorized = matmul(d,packed)
   expected = reshape(a,[9])
   call assert_close(vectorized,expected,1.0e-14_dp)

   kron = kronecker_product(a(1:2,1:2),a(1:2,1:2))
   if (any(shape(kron) /= [4,4])) error stop 1
   if (abs(kron(4,4)-16.0_dp) > 1.0e-14_dp) error stop 1

   print '(a)', 'test_linalg: PASS'

contains

   subroutine assert_close(actual, reference, tolerance)
      real(dp), intent(in) :: actual(:), reference(:), tolerance
      if (size(actual) /= size(reference)) error stop 1
      if (maxval(abs(actual-reference)) > tolerance) error stop 1
   end subroutine assert_close

end program test_linalg
