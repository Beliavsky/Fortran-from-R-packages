! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
program test_fda
   use fda
   implicit none
   integer :: failures

   failures = 0
   call test_bases(failures)
   call test_numeric(failures)
   call test_fd_smoothing(failures)
   call test_pca_cca(failures)
   call test_ode(failures)
   if (failures /= 0) then
      write (*, '(a,i0)') 'failed checks: ', failures
      error stop 1
   end if
   write (*, '(a)') 'all fda tests passed'

contains

   subroutine test_bases(failures)
      integer, intent(inout) :: failures !! Running failure count incremented for each unsuccessful basis check.
      type(basis_type) :: basis
      real(dp), allocatable :: mat(:, :), penalty(:, :)
      real(dp) :: breaks(3), expected(5, 4), exponents(3), x(5)
      integer :: info

      x = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
      breaks = [0.0_dp, 0.5_dp, 1.0_dp]
      expected = reshape([ &
         1.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.625_dp, 0.5_dp, 0.125_dp, 0.0_dp, &
         0.0_dp, 0.125_dp, 0.5_dp, 0.625_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.25_dp, 1.0_dp], [5, 4])
      call make_bspline_basis(breaks, 3, basis, info)
      call assert_int(info, 0, 'B-spline constructor', failures)
      call eval_basis(x, basis, 0, mat, info)
      call assert_int(info, 0, 'B-spline evaluation', failures)
      call assert_matrix(mat, expected, 2.0e-12_dp, 'B-spline values', failures)
      call assert_vector(sum(mat, dim=2), spread(1.0_dp, 1, size(x)), 2.0e-12_dp, &
         'B-spline partition of unity', failures)
      call eval_basis(x, basis, 1, mat, info)
      call assert_vector(sum(mat, dim=2), spread(0.0_dp, 1, size(x)), 5.0e-12_dp, &
         'B-spline derivative partition', failures)

      call make_fourier_basis([0.0_dp, 1.0_dp], 3, 1.0_dp, basis, info)
      call assert_int(info, 0, 'Fourier constructor', failures)
      call eval_basis([0.0_dp], basis, 0, mat, info)
      call assert_vector(mat(1, :), [1.0_dp, 0.0_dp, sqrt(2.0_dp)], 2.0e-12_dp, &
         'Fourier value at zero', failures)
      call eval_basis([0.0_dp], basis, 1, mat, info)
      call assert_close(mat(1, 2), sqrt(2.0_dp) * 2.0_dp * acos(-1.0_dp), 2.0e-11_dp, &
         'Fourier first derivative', failures)
      call basis_penalty(basis, 0, penalty, info, 65)
      call assert_matrix(penalty, identity_matrix(3), 2.0e-6_dp, 'Fourier Gram matrix', failures)

      exponents = [0.0_dp, 1.0_dp, 2.0_dp]
      call make_monomial_basis([0.0_dp, 1.0_dp], exponents, [0.0_dp, 1.0_dp], basis, info)
      call eval_basis([0.5_dp], basis, 2, mat, info)
      call assert_vector(mat(1, :), [0.0_dp, 0.0_dp, 2.0_dp], 2.0e-12_dp, &
         'monomial second derivative', failures)
      call basis_penalty(basis, 0, penalty, info, 65)
      expected = 0.0_dp
      call assert_close(penalty(1, 1), 1.0_dp, 2.0e-9_dp, 'monomial Gram 11', failures)
      call assert_close(penalty(1, 2), 0.5_dp, 2.0e-9_dp, 'monomial Gram 12', failures)
      call assert_close(penalty(2, 2), 1.0_dp / 3.0_dp, 2.0e-9_dp, 'monomial Gram 22', failures)
      call basis_penalty(basis, 2, penalty, info, 65)
      call assert_close(penalty(3, 3), 4.0_dp, 2.0e-9_dp, 'monomial D2 penalty', failures)

      call make_exponential_basis([0.0_dp, 1.0_dp], [2.0_dp], basis, info)
      call eval_basis([0.25_dp], basis, 2, mat, info)
      call assert_close(mat(1, 1), 4.0_dp * exp(0.5_dp), 2.0e-12_dp, 'exponential derivative', failures)

      call make_power_basis([0.25_dp, 1.0_dp], [0.5_dp, 2.5_dp], basis, info)
      call eval_basis([1.0_dp], basis, 1, mat, info)
      call assert_vector(mat(1, :), [0.5_dp, 2.5_dp], 2.0e-12_dp, 'power derivative', failures)

      call make_polygonal_basis(breaks, basis, info)
      call eval_basis([0.25_dp], basis, 0, mat, info)
      call assert_vector(mat(1, :), [0.5_dp, 0.5_dp, 0.0_dp], 2.0e-12_dp, 'polygonal basis', failures)
   end subroutine test_bases

   subroutine test_numeric(failures)
      integer, intent(inout) :: failures !! Running failure count incremented for each unsuccessful numerical-helper check.
      real(dp), allocatable :: dy(:), lmat(:, :), mmat(:, :), points(:), result(:, :), values(:), weights(:), yint(:)
      real(dp) :: amat(2, 2), bmat(2, 2), cmat(2, 2), xa(3), xmat(3, 1), ya(3, 1), ymat(3, 1)
      integer :: info

      xmat(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp]
      ymat(:, 1) = [1.0_dp, 1.0_dp, 1.0_dp]
      call trapz_mat(xmat, ymat, 1.0_dp, result=result, info=info)
      call assert_close(result(1, 1), 2.0_dp, 2.0e-12_dp, 'trapz matrix', failures)
      call quadset([0.0_dp, 0.5_dp, 1.0_dp], 5, points, weights, info)
      call assert_close(sum(weights), 1.0_dp, 2.0e-12_dp, 'Simpson weight sum', failures)

      xa = [0.0_dp, 1.0_dp, 2.0_dp]
      ya(:, 1) = xa * xa
      call polint_matrix(xa, ya, 1.5_dp, yint, dy, info)
      call assert_close(yint(1), 2.25_dp, 2.0e-12_dp, 'polynomial interpolation', failures)

      amat = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
      bmat = reshape([1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      call symsolve(amat, bmat, result, info)
      call assert_matrix(matmul(amat, result), bmat, 2.0e-11_dp, 'symmetric solve', failures)

      amat = 0.0_dp
      bmat = 0.0_dp
      cmat = 0.0_dp
      amat(1, 1) = 2.0_dp
      amat(2, 2) = 1.0_dp
      bmat(1, 1) = 1.0_dp
      bmat(2, 2) = 1.0_dp
      cmat = bmat
      call geigen(amat, bmat, cmat, values, lmat, mmat, info)
      call assert_vector(values, [2.0_dp, 1.0_dp], 2.0e-10_dp, 'generalized singular values', failures)
      if (.not. zero_find([-2.0_dp, 3.0_dp])) then
         failures = failures + 1
         write (*, '(a)') 'FAIL zero_find bracket'
      end if
      if (zero_find([2.0_dp, 3.0_dp])) then
         failures = failures + 1
         write (*, '(a)') 'FAIL zero_find no bracket'
      end if
   end subroutine test_numeric

   subroutine test_fd_smoothing(failures)
      integer, intent(inout) :: failures !! Running failure count incremented for functional-data and smoothing checks.
      type(basis_type) :: basis
      type(fd_type) :: fdobj
      type(smooth_result_type) :: fit
      real(dp), allocatable :: inner(:, :), values(:, :)
      real(dp) :: coefs(3, 2), df, lambda, t(11), y(11, 1)
      integer :: i, info

      call make_monomial_basis([0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp, 2.0_dp], [0.0_dp, 1.0_dp], basis, info)
      coefs(:, 1) = [1.0_dp, 2.0_dp, 1.0_dp]
      coefs(:, 2) = [0.0_dp, -1.0_dp, 2.0_dp]
      call make_fd(coefs, basis, fdobj, info)
      call eval_fd([0.0_dp, 0.5_dp, 1.0_dp], fdobj, 0, values, info)
      call assert_vector(values(:, 1), [1.0_dp, 2.25_dp, 4.0_dp], 2.0e-12_dp, 'fd evaluation', failures)
      call inprod_fd(fdobj, fdobj, 0, 0, inner, info, 65)
      if (info /= 0 .or. maxval(abs(inner - transpose(inner))) > 2.0e-10_dp) then
         failures = failures + 1
         write (*, '(a)') 'FAIL fd inner product symmetry'
      end if

      do i = 1, size(t)
         t(i) = real(i - 1, dp) / real(size(t) - 1, dp)
         y(i, 1) = 1.0_dp + 2.0_dp * t(i) + t(i) * t(i)
      end do
      call smooth_basis(t, y, basis, 0.0_dp, 2, fit, info)
      call assert_int(info, 0, 'unpenalized smoothing', failures)
      call assert_vector(fit%fd%coefs(:, 1), [1.0_dp, 2.0_dp, 1.0_dp], 2.0e-9_dp, &
         'exact polynomial smoothing coefficients', failures)
      call assert_close(fit%df, 3.0_dp, 2.0e-9_dp, 'unpenalized smoothing df', failures)
      call assert_close(fit%sse, 0.0_dp, 2.0e-16_dp, 'exact polynomial smoothing SSE', failures)

      call lambda_to_df(t, basis, 2, 0.1_dp, df, info)
      if (info /= 0 .or. df >= 3.0_dp .or. df <= 2.0_dp) then
         failures = failures + 1
         write (*, '(a,f12.6)') 'FAIL penalized df: ', df
      end if
      call df_to_lambda(t, basis, 2, df, lambda, info)
      call assert_close(lambda, 0.1_dp, 2.0e-6_dp, 'df-to-lambda inversion', failures)
   end subroutine test_fd_smoothing

   subroutine test_pca_cca(failures)
      integer, intent(inout) :: failures !! Running failure count incremented for functional PCA and CCA checks.
      type(basis_type) :: basis
      type(cca_result_type) :: cca
      type(fd_type) :: fdobj
      type(pca_result_type) :: pca
      real(dp) :: coefs(3, 4)
      integer :: info

      call make_fourier_basis([0.0_dp, 1.0_dp], 3, 1.0_dp, basis, info)
      coefs = 0.0_dp
      coefs(2, :) = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
      coefs(3, :) = [0.5_dp, 0.5_dp, -0.5_dp, -0.5_dp]
      call make_fd(coefs, basis, fdobj, info)
      call pca_fd(fdobj, 2, basis, 0.0_dp, 2, pca, info)
      call assert_int(info, 0, 'functional PCA', failures)
      call assert_close(pca%varprop(1), 0.8_dp, 2.0e-5_dp, 'PCA first variance proportion', failures)
      call assert_close(pca%varprop(2), 0.2_dp, 2.0e-5_dp, 'PCA second variance proportion', failures)

      call cca_fd(fdobj, fdobj, 2, 1.0e-8_dp, 1.0e-8_dp, 0, 0, cca, info)
      call assert_int(info, 0, 'functional CCA', failures)
      if (cca%correlations(1) < 0.999_dp .or. cca%correlations(2) < 0.999_dp) then
         failures = failures + 1
         write (*, '(a,2f12.6)') 'FAIL CCA correlations: ', cca%correlations(1:2)
      end if
   end subroutine test_pca_cca

   subroutine test_ode(failures)
      integer, intent(inout) :: failures !! Running failure count incremented for adaptive linear-ODE solver checks.
      type(basis_type) :: basis
      type(fd_type) :: weights(2)
      type(ode_solution_type) :: solution
      real(dp) :: coef0(1, 1), coef1(1, 1), ystart(2, 1)
      integer :: info, last

      call make_constant_basis([0.0_dp, 0.5_dp * acos(-1.0_dp)], basis, info)
      coef1(1, 1) = 1.0_dp
      coef0(1, 1) = 0.0_dp
      call make_fd(coef1, basis, weights(1), info)
      call make_fd(coef0, basis, weights(2), info)
      ystart(:, 1) = [1.0_dp, 0.0_dp]
      call odesolv(weights, ystart, solution, info, eps=1.0e-8_dp)
      call assert_int(info, 0, 'adaptive ODE solve', failures)
      last = size(solution%t)
      call assert_close(solution%t(last), 0.5_dp * acos(-1.0_dp), 2.0e-12_dp, 'ODE final time', failures)
      call assert_close(solution%y(1, 1, last), 0.0_dp, 2.0e-7_dp, 'ODE cosine solution', failures)
      call assert_close(solution%y(2, 1, last), -1.0_dp, 2.0e-7_dp, 'ODE derivative solution', failures)
   end subroutine test_ode


   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Positive order of the identity matrix to construct.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine assert_close(actual, expected, tol, label, failures)
      real(dp), intent(in) :: actual !! Scalar value produced by the routine under test.
      real(dp), intent(in) :: expected !! Independent reference value used for comparison.
      real(dp), intent(in) :: tol !! Maximum absolute error accepted by this deterministic check.
      character(len=*), intent(in) :: label !! Human-readable check name printed when the comparison fails.
      integer, intent(inout) :: failures !! Running failure count incremented when the comparison exceeds tolerance.

      if (abs(actual - expected) > tol) then
         failures = failures + 1
         write (*, '(a,2es20.10)') 'FAIL '//trim(label)//': ', actual, expected
      end if
   end subroutine assert_close

   subroutine assert_int(actual, expected, label, failures)
      integer, intent(in) :: actual !! Integer status or count produced by the routine under test.
      integer, intent(in) :: expected !! Expected integer status or count.
      character(len=*), intent(in) :: label !! Human-readable check name printed on mismatch.
      integer, intent(inout) :: failures !! Running failure count incremented when the integers differ.

      if (actual /= expected) then
         failures = failures + 1
         write (*, '(a,2i8)') 'FAIL '//trim(label)//': ', actual, expected
      end if
   end subroutine assert_int

   subroutine assert_vector(actual, expected, tol, label, failures)
      real(dp), intent(in) :: actual(:) !! Vector produced by the routine under test.
      real(dp), intent(in) :: expected(:) !! Independent reference vector with matching shape.
      real(dp), intent(in) :: tol !! Maximum componentwise absolute error accepted.
      character(len=*), intent(in) :: label !! Human-readable check name printed on mismatch.
      integer, intent(inout) :: failures !! Running failure count incremented when shape or values differ.

      if (size(actual) /= size(expected)) then
         failures = failures + 1
         write (*, '(a)') 'FAIL '//trim(label)//' shape'
      else if (maxval(abs(actual - expected)) > tol) then
         failures = failures + 1
         write (*, '(a,es20.10)') 'FAIL '//trim(label)//' max error: ', maxval(abs(actual - expected))
      end if
   end subroutine assert_vector

   subroutine assert_matrix(actual, expected, tol, label, failures)
      real(dp), intent(in) :: actual(:, :) !! Matrix produced by the routine under test.
      real(dp), intent(in) :: expected(:, :) !! Independent reference matrix with matching shape.
      real(dp), intent(in) :: tol !! Maximum elementwise absolute error accepted.
      character(len=*), intent(in) :: label !! Human-readable check name printed on mismatch.
      integer, intent(inout) :: failures !! Running failure count incremented when shape or values differ.

      if (any(shape(actual) /= shape(expected))) then
         failures = failures + 1
         write (*, '(a)') 'FAIL '//trim(label)//' shape'
      else if (maxval(abs(actual - expected)) > tol) then
         failures = failures + 1
         write (*, '(a,es20.10)') 'FAIL '//trim(label)//' max error: ', maxval(abs(actual - expected))
      end if
   end subroutine assert_matrix

end program test_fda
