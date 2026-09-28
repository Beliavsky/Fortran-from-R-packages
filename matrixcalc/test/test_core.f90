program test_core
   use matrixcalc, only : dp, direct_prod, direct_sum, hadamard_prod, frobenius_prod, &
      set_submatrix, shift_down, shift_up, shift_left, shift_right, vec, vech, &
      lower_triangle, upper_triangle, matrix_trace, entrywise_norm, frobenius_norm, &
      hilbert_schmidt_norm, maximum_norm, one_norm, inf_norm
   implicit none

   real(dp) :: a(2, 2)
   real(dp) :: b(2, 2)
   real(dp), allocatable :: z(:, :)
   real(dp), allocatable :: v(:)
   real(dp) :: expected4(4, 4)
   real(dp) :: expected3(3)
   integer :: i

   a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [2, 2])
   b = reshape([0.0_dp, 6.0_dp, 5.0_dp, 7.0_dp], [2, 2])

   z = direct_prod(a, b)
   expected4 = reshape([ &
      0.0_dp, 6.0_dp, 0.0_dp, 18.0_dp, &
      5.0_dp, 7.0_dp, 15.0_dp, 21.0_dp, &
      0.0_dp, 12.0_dp, 0.0_dp, 24.0_dp, &
      10.0_dp, 14.0_dp, 20.0_dp, 28.0_dp], [4, 4])
   call assert_matrix_close(z, expected4, 1.0e-12_dp, 'direct_prod')

   z = direct_sum(a, b)
   if (size(z, 1) /= 4 .or. size(z, 2) /= 4) error stop 'direct_sum shape'
   call assert_matrix_close(z(1:2, 1:2), a, 1.0e-12_dp, 'direct_sum A block')
   call assert_matrix_close(z(3:4, 3:4), b, 1.0e-12_dp, 'direct_sum B block')
   if (maxval(abs(z(1:2, 3:4))) > 1.0e-12_dp) error stop 'direct_sum upper off block'
   if (maxval(abs(z(3:4, 1:2))) > 1.0e-12_dp) error stop 'direct_sum lower off block'

   z = hadamard_prod(a, b)
   call assert_matrix_close(z, a * b, 1.0e-12_dp, 'hadamard_prod')
   call assert_close(frobenius_prod(a, b), sum(a * b), 1.0e-12_dp, 'frobenius_prod')

   z = set_submatrix(reshape([(real(i, dp), i = 1, 16)], [4, 4]), &
      reshape([100.0_dp, 101.0_dp, 102.0_dp, 103.0_dp], [2, 2]), 2, 2)
   call assert_close(z(2, 2), 100.0_dp, 1.0e-12_dp, 'set_submatrix 11')
   call assert_close(z(3, 3), 103.0_dp, 1.0e-12_dp, 'set_submatrix 22')

   z = shift_down(a, 1, -1.0_dp)
   call assert_matrix_close(z, reshape([-1.0_dp, 1.0_dp, -1.0_dp, 2.0_dp], [2, 2]), &
      1.0e-12_dp, 'shift_down')
   z = shift_up(a, 1, -1.0_dp)
   call assert_matrix_close(z, reshape([3.0_dp, -1.0_dp, 4.0_dp, -1.0_dp], [2, 2]), &
      1.0e-12_dp, 'shift_up')
   z = shift_left(a, 1, -1.0_dp)
   call assert_matrix_close(z, reshape([2.0_dp, 4.0_dp, -1.0_dp, -1.0_dp], [2, 2]), &
      1.0e-12_dp, 'shift_left')
   z = shift_right(a, 1, -1.0_dp)
   call assert_matrix_close(z, reshape([-1.0_dp, -1.0_dp, 1.0_dp, 3.0_dp], [2, 2]), &
      1.0e-12_dp, 'shift_right')

   v = vec(a)
   call assert_vector_close(v, [1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], 1.0e-12_dp, 'vec')
   v = vech(reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, &
      3.0_dp, 5.0_dp, 6.0_dp], [3, 3]))
   call assert_vector_close(v, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], &
      1.0e-12_dp, 'vech')

   z = lower_triangle(a)
   call assert_matrix_close(z, reshape([1.0_dp, 3.0_dp, 0.0_dp, 4.0_dp], [2, 2]), &
      1.0e-12_dp, 'lower_triangle')
   z = upper_triangle(a)
   call assert_matrix_close(z, reshape([1.0_dp, 0.0_dp, 2.0_dp, 4.0_dp], [2, 2]), &
      1.0e-12_dp, 'upper_triangle')
   call assert_close(matrix_trace(a), 5.0_dp, 1.0e-12_dp, 'matrix_trace')

   expected3 = [3.0_dp, -4.0_dp, 12.0_dp]
   call assert_close(entrywise_norm(expected3, 2.0_dp), 13.0_dp, 1.0e-12_dp, 'entrywise_norm')
   call assert_close(frobenius_norm(a), sqrt(30.0_dp), 1.0e-12_dp, 'frobenius_norm')
   call assert_close(hilbert_schmidt_norm(a), sqrt(30.0_dp), 1.0e-12_dp, 'hilbert_schmidt_norm')
   call assert_close(maximum_norm(a), 4.0_dp, 1.0e-12_dp, 'maximum_norm')
   call assert_close(one_norm(a), 6.0_dp, 1.0e-12_dp, 'one_norm')
   call assert_close(inf_norm(a), 7.0_dp, 1.0e-12_dp, 'inf_norm')

   print '(a)', 'test_core: PASS'

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

   subroutine assert_vector_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual(:) !! Computed vector.
      real(dp), intent(in) :: expected(:) !! Reference vector.
      real(dp), intent(in) :: tol !! Allowed maximum absolute error.
      character(len=*), intent(in) :: label !! Test label reported on failure.

      if (size(actual) /= size(expected)) error stop 'vector assertion shape failed'
      if (maxval(abs(actual - expected)) > tol) then
         print '(a,1x,es24.16)', trim(label), maxval(abs(actual - expected))
         error stop 'vector assertion failed'
      end if
   end subroutine assert_vector_close

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

end program test_core
