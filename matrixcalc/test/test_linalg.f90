program test_linalg
   use matrixcalc, only : dp, lu_result, lu_decomposition, matrix_inverse, matrix_power, &
      matrix_rank, spectral_norm, svd_inverse, is_square_matrix, is_symmetric_matrix, &
      is_skew_symmetric_matrix, is_diagonal_matrix, is_idempotent_matrix, &
      is_singular_matrix, is_non_singular_matrix, is_positive_definite, &
      is_positive_semi_definite, is_negative_definite, is_negative_semi_definite, &
      is_indefinite
   implicit none

   real(dp) :: a(3, 3)
   real(dp) :: diag3(3, 3)
   real(dp) :: skew(2, 2)
   real(dp), allocatable :: inv(:, :)
   real(dp), allocatable :: pwr(:, :)
   real(dp), allocatable :: pinv(:, :)
   type(lu_result) :: lu
   integer :: status

   a = reshape([4.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 5.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp, 3.0_dp], [3, 3])
   diag3 = 0.0_dp
   diag3(1, 1) = 1.0_dp
   diag3(2, 2) = 2.0_dp
   diag3(3, 3) = 3.0_dp
   skew = reshape([0.0_dp, -2.0_dp, 2.0_dp, 0.0_dp], [2, 2])

   if (.not. is_square_matrix(a)) error stop 'is_square_matrix'
   if (.not. is_symmetric_matrix(a)) error stop 'is_symmetric_matrix'
   if (.not. is_skew_symmetric_matrix(skew)) error stop 'is_skew_symmetric_matrix'
   if (.not. is_diagonal_matrix(diag3)) error stop 'is_diagonal_matrix'

   if (.not. is_idempotent_matrix(reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]))) then
      error stop 'is_idempotent_matrix'
   end if
   if (.not. is_non_singular_matrix(a, status=status)) error stop 'is_non_singular_matrix'
   if (status /= 0) error stop 'is_non_singular status'
   if (.not. is_singular_matrix(reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2]))) then
      error stop 'is_singular_matrix'
   end if

   if (.not. is_positive_definite(a, status=status)) error stop 'is_positive_definite'
   if (status /= 0) error stop 'positive definite status'
   if (.not. is_positive_semi_definite(reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]))) then
      error stop 'is_positive_semi_definite'
   end if
   if (.not. is_negative_definite(-a)) error stop 'is_negative_definite'
   if (.not. is_negative_semi_definite(-reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]))) then
      error stop 'is_negative_semi_definite'
   end if
   if (.not. is_indefinite(reshape([1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2]))) then
      error stop 'is_indefinite'
   end if

   lu = lu_decomposition(a, status)
   if (status /= 0) error stop 'lu_decomposition status'
   call assert_matrix_close(matmul(lu%l, lu%u), a, 1.0e-11_dp, 'LU reconstruction')

   inv = matrix_inverse(a, status)
   if (status /= 0) error stop 'matrix_inverse status'
   call assert_matrix_close(matmul(a, inv), identity_matrix(3), 1.0e-11_dp, 'inverse identity')

   pwr = matrix_power(a, 2, status)
   if (status /= 0) error stop 'matrix_power positive status'
   call assert_matrix_close(pwr, matmul(a, a), 1.0e-11_dp, 'matrix_power positive')
   pwr = matrix_power(a, -1, status)
   if (status /= 0) error stop 'matrix_power negative status'
   call assert_matrix_close(pwr, inv, 1.0e-11_dp, 'matrix_power negative')

   if (matrix_rank(diag3, status) /= 3 .or. status /= 0) error stop 'matrix_rank full'
   if (matrix_rank(reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2]), status) /= 1) then
      error stop 'matrix_rank deficient'
   end if
   if (status /= 0) error stop 'matrix_rank deficient status'

   call assert_close(spectral_norm(diag3, status), 3.0_dp, 1.0e-11_dp, 'spectral_norm')
   if (status /= 0) error stop 'spectral_norm status'
   call assert_close(spectral_norm([3.0_dp, 4.0_dp]), 5.0_dp, 1.0e-12_dp, 'spectral_norm vector')

   pinv = svd_inverse(a, status)
   if (status /= 0) error stop 'svd_inverse status'
   call assert_matrix_close(pinv, inv, 1.0e-10_dp, 'svd_inverse')

   print '(a)', 'test_linalg: PASS'

contains

   pure function identity_matrix(n) result(aout)
      integer, intent(in) :: n !! Order of the identity matrix.
      real(dp), allocatable :: aout(:, :)
      integer :: i

      allocate(aout(n, n))
      aout = 0.0_dp
      do i = 1, n
         aout(i, i) = 1.0_dp
      end do
   end function identity_matrix

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

end program test_linalg
