program test_multi_rhs
    use rlinsolve, only : dp, solver_matrix_result, sparse_matrix, sparse_from_dense, &
        lsolve_jacobi, lsolve_gs, lsolve_sor, lsolve_ssor, lsolve_cg, lsolve_bicg, &
        lsolve_bicgstab, lsolve_cgs, lsolve_cheby, lsolve_gmres, lsolve_qmr
    implicit none

    real(dp) :: a(3, 3)
    real(dp) :: b(3, 2)
    real(dp) :: x0(3)
    real(dp) :: xtrue(3, 2)
    type(sparse_matrix) :: as
    type(solver_matrix_result) :: out

    a = reshape([4.0_dp, 1.0_dp, 0.0_dp, &
                 1.0_dp, 3.0_dp, 1.0_dp, &
                 0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
    xtrue(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
    xtrue(:, 2) = [-1.0_dp, 1.0_dp, 0.5_dp]
    b = matmul(a, xtrue)
    x0 = 0.0_dp
    as = sparse_from_dense(a)

    call lsolve_jacobi(a, b, out, xinit=x0, reltol=1.0e-9_dp, maxiter=5000)
    call check_many(out, xtrue, 1.0e-6_dp, 'multi jacobi')

    call lsolve_gs(as, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=2000)
    call check_many(out, xtrue, 1.0e-7_dp, 'multi gs')

    call lsolve_sor(as, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=2000, weight=1.0_dp)
    call check_many(out, xtrue, 1.0e-6_dp, 'multi sor')

    call lsolve_ssor(a, b, out, xinit=x0, reltol=1.0e-9_dp, maxiter=2000, weight=1.0_dp)
    call check_many(out, xtrue, 1.0e-6_dp, 'multi ssor')

    call lsolve_cg(a, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check_many(out, xtrue, 1.0e-7_dp, 'multi cg')

    call lsolve_bicg(as, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check_many(out, xtrue, 2.0e-6_dp, 'multi bicg')

    call lsolve_bicgstab(a, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check_many(out, xtrue, 2.0e-6_dp, 'multi bicgstab')

    call lsolve_cgs(as, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check_many(out, xtrue, 2.0e-6_dp, 'multi cgs')

    call lsolve_cheby(a, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=500, &
        eigmin=1.2679491924311228_dp, eigmax=4.732050807568877_dp)
    call check_many(out, xtrue, 2.0e-6_dp, 'multi cheby')

    call lsolve_gmres(as, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=100, restart=2)
    call check_many(out, xtrue, 1.0e-6_dp, 'multi gmres')

    call lsolve_qmr(a, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=100)
    call check_many(out, xtrue, 2.0e-6_dp, 'multi qmr')

    print '(a)', 'test_multi_rhs: PASS'

contains

    subroutine check_many(result, expected, tolerance, name)
        type(solver_matrix_result), intent(in) :: result !! Multi-column result to check for convergence and accuracy.
        real(dp), intent(in) :: expected(:, :) !! Expected solution matrix with one right-hand side per column.
        real(dp), intent(in) :: tolerance !! Maximum accepted absolute solution error.
        character(len=*), intent(in) :: name !! Method label included in a failing diagnostic.

        if (.not. all(result%converged)) then
            print *, trim(name), ' flags=', result%flag, ' iter=', result%iter
            error stop 'multi-RHS convergence'
        end if
        if (maxval(abs(result%x - expected)) > tolerance) then
            print *, trim(name), result%x
            error stop 'multi-RHS solution'
        end if
        if (size(result%column) /= size(expected, 2)) error stop 'multi-RHS histories'
    end subroutine check_many

end program test_multi_rhs
