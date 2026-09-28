program test_iterative
   use sanic, only: dp, sparse_matrix, sparsify, solve_cg_sparse, sanic_cg, sanic_bicgstab, sanic_lscg
   implicit none

   real(dp), parameter :: tol = 1.0e-8_dp
   real(dp) :: aspd(3, 3), agen(3, 3), arect(4, 3), b3(3, 1), b4(4, 1), truth(3, 1)
   real(dp) :: b_multi(3, 2), truth_multi(3, 2)
   real(dp) :: a5(5, 5), b5(5, 1), truth5(5, 1)
   real(dp), allocatable :: x(:, :)
   type(sparse_matrix) :: a
   integer :: info, iterations
   real(dp) :: err
   logical :: converged

   truth(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp]
   aspd = reshape([4.0_dp, 1.0_dp, 0.0_dp, &
                   1.0_dp, 3.0_dp, 1.0_dp, &
                   0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
   b3 = matmul(aspd, truth)
   a = sparsify(aspd)
   call solve_cg_sparse(a, x, info, b=b3, method=sanic_cg, iter=50, tol=1.0e-12_dp, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth)) < tol, "CG")

   call solve_cg_sparse(a, x, info, b=b3, method=sanic_cg, iter=50, tol=1.0e-12_dp, precond=2, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth)) < tol, "CG incomplete Cholesky")
   call require(iterations <= 2, "CG incomplete Cholesky convergence")

   agen = reshape([4.0_dp, -1.0_dp, 0.5_dp, &
                   2.0_dp, 5.0_dp, 1.0_dp, &
                   0.0_dp, 1.0_dp, 3.0_dp], [3, 3])
   b3 = matmul(agen, truth)
   a = sparsify(agen)
   call solve_cg_sparse(a, x, info, b=b3, method=sanic_bicgstab, iter=100, tol=1.0e-12_dp, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth)) < 2.0e-8_dp, "BiCGSTAB")

   call solve_cg_sparse(a, x, info, b=b3, method=sanic_bicgstab, iter=1, tol=1.0e-30_dp, precond=0, &
                        iterations=iterations, error=err, converged=converged)
   call require(iterations == 1 .and. .not. converged .and. info == 2, "explicit iteration budget honored")

   truth_multi(:, 1) = truth(:, 1)
   truth_multi(:, 2) = [-1.0_dp, 0.25_dp, 2.0_dp]
   b_multi = matmul(agen, truth_multi)
   call solve_cg_sparse(a, x, info, b=b_multi, method=sanic_bicgstab, iter=100, tol=1.0e-12_dp, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. all(shape(x) == [3, 2]), "multiple RHS shape")
   call require(maxval(abs(x - truth_multi)) < 2.0e-8_dp, "multiple RHS values")

   call solve_cg_sparse(a, x, info, b=b3, method=sanic_bicgstab, iter=-1, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 1 .and. size(x) == 0, "negative iteration count rejected")
   call solve_cg_sparse(a, x, info, b=b3, method=sanic_bicgstab, tol=-1.0_dp, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 1 .and. size(x) == 0, "negative tolerance rejected")

   agen = reshape([4.0_dp, -1.0_dp, 0.0_dp, &
                   1.0_dp, 4.0_dp, -1.0_dp, &
                   0.0_dp, 1.0_dp, 3.0_dp], [3, 3])
   b3 = matmul(agen, truth)
   a = sparsify(agen)
   call solve_cg_sparse(a, x, info, b=b3, method=sanic_bicgstab, iter=50, tol=1.0e-12_dp, precond=2, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth)) < 2.0e-8_dp, "BiCGSTAB ILUT")
   call require(iterations <= 2, "BiCGSTAB ILUT convergence")


   truth5(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -1.0_dp]
   a5 = 0.0_dp
   do iterations = 1, 5
      a5(iterations, iterations) = 5.0_dp
   end do
   do iterations = 1, 4
      a5(iterations, iterations + 1) = 1.0_dp
      a5(iterations + 1, iterations) = 1.0_dp
   end do
   a5(1, 5) = 1.0_dp
   a5(5, 1) = 1.0_dp
   b5 = matmul(a5, truth5)
   a = sparsify(a5)
   call solve_cg_sparse(a, x, info, b=b5, method=sanic_cg, iter=50, tol=1.0e-12_dp, precond=2, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth5)) < tol, &
                "CG limited-memory incomplete Cholesky")
   call require(iterations <= 4, "CG limited-memory incomplete Cholesky convergence")

   a5 = 0.0_dp
   do iterations = 1, 5
      a5(iterations, iterations) = 5.0_dp
   end do
   a5(1, 2) = 2.0_dp
   a5(2, 3) = -1.0_dp
   a5(3, 4) = 1.5_dp
   a5(4, 5) = -0.75_dp
   a5(5, 1) = 1.25_dp
   a5(2, 1) = -0.5_dp
   a5(3, 1) = 0.3_dp
   a5(4, 2) = 0.4_dp
   a5(5, 3) = -0.2_dp
   b5 = matmul(a5, truth5)
   a = sparsify(a5)
   call solve_cg_sparse(a, x, info, b=b5, method=sanic_bicgstab, iter=50, tol=1.0e-12_dp, precond=2, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth5)) < 2.0e-8_dp, &
                "BiCGSTAB dual-threshold ILUT")
   call require(iterations <= 4, "BiCGSTAB dual-threshold ILUT convergence")

   arect = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                    0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
                    2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp], [4, 3])
   b4 = matmul(arect, truth)
   a = sparsify(arect)
   call solve_cg_sparse(a, x, info, b=b4, method=sanic_lscg, iter=100, tol=1.0e-12_dp, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth)) < 2.0e-8_dp, "LSCG")

   call solve_cg_sparse(a, x, info, b=b4, method=sanic_lscg, iter=100, tol=1.0e-12_dp, precond=2, &
                        iterations=iterations, error=err, converged=converged)
   call require(info == 0 .and. converged .and. maxval(abs(x - truth)) < 2.0e-8_dp, "LSCG precond=2 default diagonal")

   print '(a)', 'test_iterative: PASS'

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition !! True when the regression condition has passed.
      character(len=*), intent(in) :: label !! Short diagnostic label printed on failure.
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAILED:', trim(label)
         error stop 1
      end if
   end subroutine require

end program test_iterative
