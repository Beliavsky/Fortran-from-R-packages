program test_optr
   use optr, only : dp, optr_success, optr_invalid_shape, cholesky_decompose, &
      conjugate_gradient, finite_difference_jacobian, gauss_seidel, gaussian_solve, &
      hat_matrix, invert_matrix, lu_decompose, lu_solve, lu_split, newton_raphson, &
      optr_fit_design, optr_solve, predict_linear, vector_function
   implicit none

   call test_direct_solvers()
   call test_cholesky()
   call test_iterative_solvers()
   call test_inverse_lu_and_hat()
   call test_model_fit_and_predict()
   call test_jacobian_and_newton()
   call test_validation()

   print '(a)', 'test_optr: PASS'

contains

   subroutine test_direct_solvers()
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: xm(:, :)
      real(dp) :: a(3, 3)
      real(dp) :: b(3)
      real(dp) :: bm(3, 2)
      integer :: status

      a = reshape([6.0_dp, -4.0_dp, 1.0_dp, -4.0_dp, 6.0_dp, -4.0_dp, &
         1.0_dp, -4.0_dp, 6.0_dp], [3, 3], order=[2, 1])
      b = [-14.0_dp, 36.0_dp, 6.0_dp]

      call gaussian_solve(a, b, x, status)
      call assert_true(status == optr_success, 'gaussian status')
      call assert_close_vector(matmul(a, x), b, 1.0e-11_dp, 'gaussian residual')

      call lu_solve(a, b, x, status)
      call assert_true(status == optr_success, 'lu status')
      call assert_close_vector(matmul(a, x), b, 1.0e-11_dp, 'lu residual')

      bm(:, 1) = b
      bm(:, 2) = [1.0_dp, 2.0_dp, 3.0_dp]
      call optr_solve(a, bm, 'LU', xm, status)
      call assert_true(status == optr_success, 'optr matrix status')
      call assert_close_matrix(matmul(a, xm), bm, 1.0e-11_dp, 'optr multi-response residual')
   end subroutine test_direct_solvers

   subroutine test_cholesky()
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: x(:)
      real(dp) :: a(3, 3)
      real(dp) :: b(3)
      integer :: status

      a = reshape([6.0_dp, -4.0_dp, 1.0_dp, -4.0_dp, 6.0_dp, -4.0_dp, &
         1.0_dp, -4.0_dp, 6.0_dp], [3, 3], order=[2, 1])
      b = [-14.0_dp, 36.0_dp, 6.0_dp]
      call cholesky_decompose(a, l, status)
      call assert_true(status == optr_success, 'cholesky factor status')
      call assert_close_matrix(matmul(l, transpose(l)), a, 1.0e-11_dp, 'cholesky reconstruction')

      call optr_solve(a, b, 'choleski', x, status)
      call assert_true(status == optr_success, 'choleski solve status')
      call assert_close_vector(matmul(a, x), b, 1.0e-11_dp, 'choleski residual')
   end subroutine test_cholesky

   subroutine test_iterative_solvers()
      real(dp), allocatable :: convergence(:)
      real(dp), allocatable :: x(:)
      real(dp) :: a(3, 3)
      real(dp) :: b(3)
      integer :: iterations
      integer :: status

      a = reshape([4.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 4.0_dp, -2.0_dp, &
         1.0_dp, -2.0_dp, 4.0_dp], [3, 3], order=[2, 1])
      b = [12.0_dp, -1.0_dp, 5.0_dp]

      call conjugate_gradient(a, b, x, status, iterations, convergence, tol=1.0e-10_dp)
      call assert_true(status == optr_success, 'cgm status')
      call assert_true(iterations > 0, 'cgm iterations')
      call assert_close_vector(matmul(a, x), b, 1.0e-8_dp, 'cgm residual')

      call gauss_seidel(a, b, x, status, iterations, tol=1.0e-10_dp)
      call assert_true(status == optr_success, 'gauss-seidel status')
      call assert_true(iterations > 0, 'gauss-seidel iterations')
      call assert_close_vector(matmul(a, x), b, 1.0e-8_dp, 'gauss-seidel residual')
   end subroutine test_iterative_solvers

   subroutine test_inverse_lu_and_hat()
      real(dp), allocatable :: a_inverse(:, :)
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: lu(:, :)
      real(dp), allocatable :: u(:, :)
      real(dp) :: a(3, 3)
      real(dp) :: identity(3, 3)
      real(dp) :: design(4, 2)
      integer, allocatable :: permutation(:)
      integer :: i
      integer :: status

      a = reshape([0.6_dp, -0.4_dp, 1.0_dp, -0.3_dp, 0.2_dp, 0.5_dp, &
         0.6_dp, -1.0_dp, 0.5_dp], [3, 3], order=[2, 1])
      identity = 0.0_dp
      do i = 1, 3
         identity(i, i) = 1.0_dp
      end do

      call invert_matrix(a, a_inverse, status)
      call assert_true(status == optr_success, 'inverse status')
      call assert_close_matrix(matmul(a, a_inverse), identity, 1.0e-10_dp, 'inverse identity')

      call lu_decompose(a, lu, permutation, status)
      call assert_true(status == optr_success, 'lu factor status')
      call lu_split(lu, l, u, status)
      call assert_true(status == optr_success, 'lu split status')
      call assert_close_matrix(matmul(l, u), a(permutation, :), 1.0e-11_dp, 'lu reconstruction')

      design = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
         1.0_dp, 3.0_dp], [4, 2], order=[2, 1])
      call hat_matrix(design, h, status)
      call assert_true(status == optr_success, 'hat status')
      call assert_close_matrix(matmul(h, h), h, 1.0e-11_dp, 'hat idempotence')
      call assert_close_matrix(h, transpose(h), 1.0e-11_dp, 'hat symmetry')
   end subroutine test_inverse_lu_and_hat

   subroutine test_model_fit_and_predict()
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: fitted(:)
      real(dp) :: design(5, 2)
      real(dp) :: response(5)
      real(dp) :: weights(5)
      integer :: status

      design = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
         1.0_dp, 3.0_dp, 1.0_dp, 4.0_dp], [5, 2], order=[2, 1])
      response = [2.0_dp, 5.0_dp, 8.0_dp, 11.0_dp, 14.0_dp]
      weights = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]

      call optr_fit_design(design, response, beta, status, method='LU')
      call assert_true(status == optr_success, 'fit status')
      call assert_close_vector(beta, [2.0_dp, 3.0_dp], 1.0e-11_dp, 'fit coefficients')

      call optr_fit_design(design, response, beta, status, method='LU', weights=weights)
      call assert_true(status == optr_success, 'weighted fit status')
      call assert_close_vector(beta, [2.0_dp, 3.0_dp], 1.0e-11_dp, 'weighted fit coefficients')

      call predict_linear(design, beta, fitted, status)
      call assert_true(status == optr_success, 'predict status')
      call assert_close_vector(fitted, response, 1.0e-11_dp, 'predictions')
   end subroutine test_model_fit_and_predict

   subroutine test_jacobian_and_newton()
      type(vector_function) :: f
      real(dp), allocatable :: f0(:)
      real(dp), allocatable :: jac(:, :)
      real(dp), allocatable :: root(:)
      real(dp) :: expected_jac(2, 2)
      real(dp) :: x(2)
      integer :: iterations
      integer :: status

      f%evaluate => nonlinear_system
      x = [1.0_dp, 2.0_dp]
      call finite_difference_jacobian(f, x, jac, f0, status)
      call assert_true(status == optr_success, 'jacobian status')
      expected_jac = reshape([2.0001_dp, 4.0001_dp, 1.0_dp, -1.0_dp], [2, 2], order=[2, 1])
      call assert_close_matrix(jac, expected_jac, 1.0e-7_dp, 'finite-difference jacobian')

      call newton_raphson(f, [1.5_dp, 0.5_dp], root, status, iterations, tol=1.0e-10_dp)
      call assert_true(status == optr_success, 'newton status')
      call assert_true(iterations > 0, 'newton iterations')
      call assert_close_vector(root, [2.0_dp, 1.0_dp], 1.0e-7_dp, 'newton root')
   end subroutine test_jacobian_and_newton

   subroutine test_validation()
      real(dp), allocatable :: x(:)
      real(dp) :: a(2, 3)
      real(dp) :: b(2)
      integer :: status

      a = 1.0_dp
      b = 1.0_dp
      call gaussian_solve(a, b, x, status)
      call assert_true(status == optr_invalid_shape, 'invalid shape detection')
   end subroutine test_validation

   pure subroutine nonlinear_system(x, fx)
      real(dp), intent(in) :: x(:) !! Two-vector at which the nonlinear test system is evaluated.
      real(dp), intent(out) :: fx(:) !! Two-vector containing the nonlinear test equations.

      fx(1) = x(1) * x(1) + x(2) * x(2) - 5.0_dp
      fx(2) = x(1) - x(2) - 1.0_dp
   end subroutine nonlinear_system

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Test condition that must evaluate true.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the assertion fails.

      if (.not. condition) then
         write (*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed vector under test.
      real(dp), intent(in) :: expected(:) !! Reference vector of the same shape.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute elementwise difference.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the assertion fails.

      if (size(actual) /= size(expected)) then
         write (*, '(a)') 'FAIL shape: ' // trim(label)
         error stop 1
      end if
      if (maxval(abs(actual - expected)) > tolerance) then
         write (*, '(a,es14.6)') 'FAIL: ' // trim(label) // ' max error=', maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_close_vector

   subroutine assert_close_matrix(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:, :) !! Computed matrix under test.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix of the same shape.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute elementwise difference.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the assertion fails.

      if (any(shape(actual) /= shape(expected))) then
         write (*, '(a)') 'FAIL shape: ' // trim(label)
         error stop 1
      end if
      if (maxval(abs(actual - expected)) > tolerance) then
         write (*, '(a,es14.6)') 'FAIL: ' // trim(label) // ' max error=', maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_close_matrix

end program test_optr
