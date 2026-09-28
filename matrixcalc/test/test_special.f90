program test_special
   use matrixcalc, only : dp, creation_matrix, fibonacci_matrix, frobenius_matrix, &
      hankel_matrix, hilbert_matrix, pascal_matrix, stirling_matrix, &
      symmetric_pascal_matrix, toeplitz_matrix, vandermonde_matrix
   implicit none

   real(dp), allocatable :: a(:, :)

   a = creation_matrix(3)
   call assert_matrix_close(a, reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [3, 3]), 1.0e-12_dp, 'creation_matrix')

   a = fibonacci_matrix(2)
   call assert_matrix_close(a, reshape([1.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, &
      2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [3, 3]), 1.0e-12_dp, 'fibonacci_matrix')

   a = frobenius_matrix(3)
   call assert_matrix_close(a, reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 1.0_dp, -3.0_dp, 3.0_dp], [3, 3]), 1.0e-12_dp, 'frobenius_matrix')

   a = hankel_matrix(3, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])
   call assert_matrix_close(a, reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 3.0_dp, &
      4.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [3, 3]), 1.0e-12_dp, 'hankel_matrix')

   a = hilbert_matrix(3)
   call assert_close(a(1, 1), 1.0_dp, 1.0e-12_dp, 'hilbert 11')
   call assert_close(a(3, 3), 0.2_dp, 1.0e-12_dp, 'hilbert 33')

   a = pascal_matrix(4)
   call assert_matrix_close(a, reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [4, 4]), 1.0e-12_dp, 'pascal_matrix')

   a = symmetric_pascal_matrix(3)
   call assert_matrix_close(a, reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
      3.0_dp, 1.0_dp, 3.0_dp, 6.0_dp], [3, 3]), 1.0e-12_dp, 'symmetric_pascal_matrix')

   a = stirling_matrix(3)
   call assert_close(a(4, 2), 1.0_dp, 1.0e-12_dp, 'stirling S31')
   call assert_close(a(4, 3), 3.0_dp, 1.0e-12_dp, 'stirling S32')
   call assert_close(a(4, 4), 1.0_dp, 1.0e-12_dp, 'stirling S33')

   a = toeplitz_matrix(3, [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp])
   call assert_matrix_close(a, reshape([30.0_dp, 40.0_dp, 50.0_dp, 20.0_dp, 30.0_dp, &
      40.0_dp, 10.0_dp, 20.0_dp, 30.0_dp], [3, 3]), 1.0e-12_dp, 'toeplitz_matrix')

   a = vandermonde_matrix([2.0_dp, 3.0_dp], 4)
   call assert_matrix_close(a, reshape([1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 9.0_dp, &
      8.0_dp, 27.0_dp], [2, 4]), 1.0e-12_dp, 'vandermonde_matrix')

   print '(a)', 'test_special: PASS'

contains

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tol !! Allowed absolute error.
      character(len=*), intent(in) :: label !! Test label reported on failure.

      if (abs(actual - expected) > tol) then
         print '(a,2(1x,es24.16))', trim(label), actual, expected
         error stop 'scalar assertion failed'
      end if
   end subroutine assert_close

   subroutine assert_matrix_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual(:, :) !! Computed matrix.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix.
      real(dp), intent(in) :: tol !! Allowed maximum absolute error.
      character(len=*), intent(in) :: label !! Test label reported on failure.

      if (any(shape(actual) /= shape(expected))) error stop 'matrix assertion shape failed'
      if (maxval(abs(actual - expected)) > tol) then
         print '(a,1x,es24.16)', trim(label), maxval(abs(actual - expected))
         error stop 'matrix assertion failed'
      end if
   end subroutine assert_matrix_close

end program test_special
