program test_structured
   use matrixcalc, only : dp, e_matrices, h_matrices, t_matrices, u_vectors_result, &
      u_vectors, k_matrix, commutation_matrix, d_matrix, duplication_matrix, &
      l_matrix, elimination_matrix, n_matrix, vec, vech
   implicit none

   real(dp) :: a(2, 3)
   real(dp) :: s(3, 3)
   real(dp), allocatable :: k(:, :)
   real(dp), allocatable :: d(:, :)
   real(dp), allocatable :: l(:, :)
   real(dp), allocatable :: nmat(:, :)
   real(dp), allocatable :: e(:, :, :, :)
   real(dp), allocatable :: h(:, :, :, :)
   real(dp), allocatable :: tmat(:, :, :, :)
   real(dp), allocatable :: v(:)
   type(u_vectors_result) :: u
   integer :: i

   a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [2, 3])
   k = k_matrix(2, 3)
   call assert_vector_close(matmul(k, vec(a)), vec(transpose(a)), 1.0e-12_dp, 'K vec identity')
   call assert_matrix_close(commutation_matrix(2, 3), k, 1.0e-12_dp, 'commutation alias')

   s = reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, &
      3.0_dp, 5.0_dp, 6.0_dp], [3, 3])
   d = d_matrix(3)
   call assert_vector_close(matmul(d, vech(s)), vec(s), 1.0e-12_dp, 'D vech identity')
   call assert_matrix_close(duplication_matrix(3), d, 1.0e-12_dp, 'duplication alias')

   l = l_matrix(3)
   call assert_vector_close(matmul(l, vec(s)), vech(s), 1.0e-12_dp, 'L vec identity')
   call assert_matrix_close(elimination_matrix(3), l, 1.0e-12_dp, 'elimination alias')

   nmat = n_matrix(3)
   v = matmul(nmat, vec(reshape([(real(i, dp), i = 1, 9)], [3, 3])))
   call assert_vector_close(v, vec(0.5_dp * (reshape([(real(i, dp), i = 1, 9)], [3, 3]) + &
      transpose(reshape([(real(i, dp), i = 1, 9)], [3, 3])))), 1.0e-12_dp, 'N symmetrizer')

   e = e_matrices(3)
   if (size(e, 1) /= 3 .or. size(e, 4) /= 3) error stop 'E shape'
   if (abs(e(2, 3, 2, 3) - 1.0_dp) > 1.0e-12_dp) error stop 'E basis entry'
   if (abs(sum(e(:, :, 2, 3)) - 1.0_dp) > 1.0e-12_dp) error stop 'E basis sum'

   h = h_matrices(2, 3)
   if (abs(h(2, 3, 2, 3) - 1.0_dp) > 1.0e-12_dp) error stop 'H basis entry'
   if (abs(sum(h(:, :, 2, 3)) - 1.0_dp) > 1.0e-12_dp) error stop 'H basis sum'

   tmat = t_matrices(3)
   if (abs(tmat(2, 3, 2, 3) - 1.0_dp) > 1.0e-12_dp) error stop 'T direct entry'
   if (abs(tmat(3, 2, 2, 3) - 1.0_dp) > 1.0e-12_dp) error stop 'T transpose entry'
   if (abs(sum(tmat(:, :, 2, 3)) - 2.0_dp) > 1.0e-12_dp) error stop 'T basis sum'

   u = u_vectors(3)
   if (any(shape(u%k) /= [3, 3])) error stop 'u k shape'
   if (u%k(1, 1) /= 1 .or. u%k(3, 1) /= 3 .or. u%k(2, 2) /= 4 .or. u%k(3, 3) /= 6) then
      error stop 'u k mapping'
   end if
   if (abs(sum(u%ident) - 6.0_dp) > 1.0e-12_dp) error stop 'u identity'

   print '(a)', 'test_structured: PASS'

contains

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

end program test_structured
