program test_rarpack
   use rarpack, only: dp, eigs, eigs_result, eigs_sym, eigs_sym_result, svds, svds_result
   implicit none

   call test_symmetric_eigs()
   call test_symmetric_shift()
   call test_general_eigs()
   call test_complex_eigenpair()
   call test_general_shift()
   call test_svd_tall()
   call test_svd_nontrivial()
   call test_svd_wide_and_vector_counts()
   call test_retvec_false()
   call test_validation()

   print '(a)', 'test_rarpack: PASS'

contains

   subroutine test_symmetric_eigs()
      real(dp) :: a(10, 10)
      type(eigs_sym_result) :: r
      integer :: i

      a = 0.0_dp
      do i = 1, 10
         a(i, i) = real(i, dp)
      end do
      a(10, 1) = 123.0_dp
      r = eigs_sym(a, 3, which = 'LA', lower = .false.)
      call require(r%info == 0, 'symmetric LA status')
      call require(r%nconv == 3, 'symmetric LA convergence count')
      call require(maxval(abs(r%values - [10.0_dp, 9.0_dp, 8.0_dp])) < 1.0e-8_dp, 'symmetric LA values')
      call require_symmetric_residual(diagonal_matrix(10), r, 1.0e-7_dp, 'symmetric LA residual')

      r = eigs_sym(diagonal_matrix(10), 3, which = 'SA')
      call require(r%info == 0, 'symmetric SA status')
      call require(maxval(abs(r%values - [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-8_dp, 'symmetric SA values')
   end subroutine test_symmetric_eigs

   subroutine test_symmetric_shift()
      real(dp) :: a(10, 10)
      type(eigs_sym_result) :: r

      a = diagonal_matrix(10)
      r = eigs_sym(a, 2, which = 'LM', sigma = 5.2_dp)
      call require(r%info == 0, 'symmetric shift status')
      call require(r%nconv == 2, 'symmetric shift convergence count')
      call require(abs(r%values(1) - 5.0_dp) < 1.0e-8_dp, 'symmetric shift nearest value')
      call require(abs(r%values(2) - 6.0_dp) < 1.0e-8_dp, 'symmetric shift second value')
      call require_symmetric_residual(a, r, 1.0e-7_dp, 'symmetric shift residual')
   end subroutine test_symmetric_shift

   subroutine test_general_eigs()
      real(dp) :: a(8, 8)
      type(eigs_result) :: r
      integer :: i

      a = 0.0_dp
      do i = 1, 8
         a(i, i) = real(i, dp)
      end do
      do i = 1, 7
         a(i, i + 1) = 0.2_dp * real(i, dp)
      end do
      r = eigs(a, 3, which = 'LR')
      call require(r%info == 0, 'general LR status')
      call require(r%nconv == 3, 'general LR convergence count')
      call require(maxval(abs(real(r%values, dp) - [8.0_dp, 7.0_dp, 6.0_dp])) < 1.0e-7_dp, 'general LR values')
      call require(maxval(abs(aimag(r%values))) < 1.0e-8_dp, 'general LR imaginary parts')
      call require_general_residual(a, r, 1.0e-6_dp, 'general LR residual')
   end subroutine test_general_eigs

   subroutine test_complex_eigenpair()
      real(dp) :: a(8, 8)
      type(eigs_result) :: r
      integer :: i

      a = 0.0_dp
      a(1, 2) = -5.0_dp
      a(2, 1) = 5.0_dp
      a(3, 4) = -3.0_dp
      a(4, 3) = 3.0_dp
      do i = 5, 8
         a(i, i) = real(i - 4, dp)
      end do

      r = eigs(a, 2, which = 'LI')
      call require(r%info == 0, 'complex eigenpair status')
      call require(r%nconv == 2, 'complex eigenpair convergence count')
      call require(abs(abs(aimag(r%values(1))) - 5.0_dp) < 1.0e-8_dp, 'complex eigenpair first value')
      call require(abs(abs(aimag(r%values(2))) - 5.0_dp) < 1.0e-8_dp, 'complex eigenpair second value')
      call require_general_residual(a, r, 1.0e-7_dp, 'complex eigenpair residual')
   end subroutine test_complex_eigenpair

   subroutine test_general_shift()
      real(dp) :: a(8, 8)
      type(eigs_result) :: r

      a = diagonal_matrix(8)
      r = eigs(a, 2, which = 'LM', sigma = 4.2_dp)
      call require(r%info == 0, 'general shift status')
      call require(r%nconv == 2, 'general shift convergence count')
      call require(abs(real(r%values(1), dp) - 4.0_dp) < 1.0e-8_dp, 'general shift nearest value')
      call require(abs(real(r%values(2), dp) - 5.0_dp) < 1.0e-8_dp, 'general shift second value')
      call require_general_residual(a, r, 1.0e-7_dp, 'general shift residual')
   end subroutine test_general_shift

   subroutine test_svd_tall()
      real(dp) :: a(6, 4)
      type(svds_result) :: r

      a = 0.0_dp
      a(1, 1) = 4.0_dp
      a(2, 2) = 3.0_dp
      a(3, 3) = 2.0_dp
      a(4, 4) = 1.0_dp
      a(5, 1) = 0.0_dp
      a(6, 2) = 0.0_dp

      r = svds(a, 2)
      call require(r%info == 0, 'tall svds status')
      call require(r%nconv == 2, 'tall svds convergence count')
      call require(maxval(abs(r%d - [4.0_dp, 3.0_dp])) < 1.0e-8_dp, 'tall svds values')
      call require_svd_residual(a, r, 1.0e-7_dp, 'tall svds residual')
   end subroutine test_svd_tall

   subroutine test_svd_nontrivial()
      real(dp) :: a(12, 7)
      real(dp), parameter :: expected(3) = [5.033286700240_dp, 3.047587050530_dp, 2.664981788334_dp]
      type(svds_result) :: r
      integer :: i
      integer :: j

      do j = 1, 7
         do i = 1, 12
            a(i, j) = sin(real(i * j, dp)) + 0.1_dp * real(i, dp) - 0.05_dp * real(j, dp)
         end do
      end do

      r = svds(a, 3)
      call require(r%info == 0, 'nontrivial svds status')
      call require(maxval(abs(r%d - expected)) < 2.0e-8_dp, 'nontrivial svds values')
      call require_svd_residual(a, r, 2.0e-7_dp, 'nontrivial svds residual')
   end subroutine test_svd_nontrivial

   subroutine test_svd_wide_and_vector_counts()
      real(dp) :: a(4, 6)
      type(svds_result) :: r

      a = 0.0_dp
      a(1, 1) = 5.0_dp
      a(2, 2) = 4.0_dp
      a(3, 3) = 2.0_dp
      a(4, 4) = 1.0_dp

      r = svds(a, 2, nu = 0, nv = 1)
      call require(r%info == 0, 'wide svds status')
      call require(maxval(abs(r%d - [5.0_dp, 4.0_dp])) < 1.0e-8_dp, 'wide svds values')
      call require(size(r%u, 2) == 0, 'wide svds nu count')
      call require(size(r%v, 2) == 1, 'wide svds nv count')
      call require(norm2(matmul(transpose(a), [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]) - 5.0_dp * r%v(:, 1)) < 1.0e-7_dp .or. &
                   norm2(matmul(transpose(a), [-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]) - 5.0_dp * r%v(:, 1)) < 1.0e-7_dp, &
                   'wide svds right vector')
   end subroutine test_svd_wide_and_vector_counts

   subroutine test_retvec_false()
      type(eigs_sym_result) :: rs
      type(eigs_result) :: rg

      rs = eigs_sym(diagonal_matrix(8), 2, retvec = .false.)
      call require(rs%info == 0, 'symmetric retvec false status')
      call require(size(rs%vectors, 2) == 0, 'symmetric retvec false vectors')

      rg = eigs(diagonal_matrix(8), 2, retvec = .false.)
      call require(rg%info == 0, 'general retvec false status')
      call require(size(rg%vectors, 2) == 0, 'general retvec false vectors')
   end subroutine test_retvec_false

   subroutine test_validation()
      type(eigs_sym_result) :: rs
      type(eigs_result) :: rg
      type(svds_result) :: rv
      real(dp) :: nonsquare(3, 4)

      nonsquare = 0.0_dp
      rg = eigs(nonsquare, 1)
      call require(rg%info < 0, 'general shape validation')
      rs = eigs_sym(diagonal_matrix(4), 4)
      call require(rs%info < 0, 'symmetric k validation')
      rv = svds(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 1, nu = 2)
      call require(rv%info < 0, 'svds vector-count validation')
   end subroutine test_validation

   function diagonal_matrix(n) result(a)
      integer, intent(in) :: n !! Order of the diagonal matrix with entries 1 through n.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = real(i, dp)
      end do
   end function diagonal_matrix

   subroutine require_symmetric_residual(a, r, tol, label)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix used to evaluate returned eigenpair residuals.
      type(eigs_sym_result), intent(in) :: r !! Symmetric eigensolver result whose vectors and values are checked.
      real(dp), intent(in) :: tol !! Maximum accepted Euclidean residual norm for each returned eigenpair.
      character(len = *), intent(in) :: label !! Diagnostic label printed if a residual exceeds tolerance.
      integer :: j
      real(dp) :: residual

      do j = 1, r%nconv
         residual = norm2(matmul(a, r%vectors(:, j)) - r%values(j) * r%vectors(:, j))
         call require(residual < tol, label)
      end do
   end subroutine require_symmetric_residual

   subroutine require_general_residual(a, r, tol, label)
      real(dp), intent(in) :: a(:, :) !! Real matrix used to evaluate returned complex eigenpair residuals.
      type(eigs_result), intent(in) :: r !! Nonsymmetric eigensolver result whose vectors and values are checked.
      real(dp), intent(in) :: tol !! Maximum accepted Euclidean residual norm for each returned eigenpair.
      character(len = *), intent(in) :: label !! Diagnostic label printed if a residual exceeds tolerance.
      integer :: j
      real(dp) :: residual

      do j = 1, r%nconv
         residual = sqrt(sum(abs(matmul(cmplx(a, 0.0_dp, kind = dp), r%vectors(:, j)) - &
                                 r%values(j) * r%vectors(:, j)) ** 2))
         call require(residual < tol, label)
      end do
   end subroutine require_general_residual

   subroutine require_svd_residual(a, r, tol, label)
      real(dp), intent(in) :: a(:, :) !! Matrix used to evaluate the returned singular triplets.
      type(svds_result), intent(in) :: r !! Truncated-SVD result whose singular vectors and values are checked.
      real(dp), intent(in) :: tol !! Maximum accepted Euclidean residual norm for each singular-vector equation.
      character(len = *), intent(in) :: label !! Diagnostic label printed if a singular residual exceeds tolerance.
      integer :: j
      real(dp) :: residual_left
      real(dp) :: residual_right

      do j = 1, min(size(r%u, 2), size(r%v, 2))
         residual_left = norm2(matmul(a, r%v(:, j)) - r%d(j) * r%u(:, j))
         residual_right = norm2(matmul(transpose(a), r%u(:, j)) - r%d(j) * r%v(:, j))
         call require(max(residual_left, residual_right) < tol, label)
      end do
   end subroutine require_svd_residual

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Assertion condition that must be true for the test to continue.
      character(len = *), intent(in) :: message !! Diagnostic description reported when the assertion fails.

      if (.not. condition) then
         write (*, '(a)') 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine require

end program test_rarpack
