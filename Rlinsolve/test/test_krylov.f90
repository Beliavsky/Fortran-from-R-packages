program test_krylov
    use rlinsolve, only : dp, sparse_matrix, solver_result, sparse_from_dense, &
        lsolve_cg, lsolve_bicg, lsolve_bicgstab, lsolve_cgs, lsolve_cheby, &
        lsolve_gmres, lsolve_qmr
    implicit none

    real(dp) :: a(3, 3)
    real(dp) :: b(3)
    real(dp) :: x0(3)
    real(dp) :: xtrue(3)
    type(sparse_matrix) :: as
    type(solver_result) :: out

    a = reshape([4.0_dp, 1.0_dp, 0.0_dp, &
                 1.0_dp, 3.0_dp, 1.0_dp, &
                 0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
    xtrue = [1.0_dp, 2.0_dp, 3.0_dp]
    b = matmul(a, xtrue)
    x0 = 0.0_dp
    as = sparse_from_dense(a)

    call lsolve_cg(a, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check(out, xtrue, 'cg')

    call lsolve_bicg(as, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check(out, xtrue, 'bicg')

    call lsolve_bicgstab(a, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check(out, xtrue, 'bicgstab')

    call lsolve_cgs(as, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100)
    call check(out, xtrue, 'cgs')

    call lsolve_cheby(a, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=500, &
        eigmin=1.2679491924311228_dp, eigmax=4.732050807568877_dp)
    call check(out, xtrue, 'cheby')

    call lsolve_cheby(a, b, out, xinit=x0, reltol=1.0e-8_dp, maxiter=500)
    call check(out, xtrue, 'cheby-estimated')

    call lsolve_gmres(as, b, out, xinit=x0, reltol=1.0e-11_dp, maxiter=100, restart=2)
    call check(out, xtrue, 'gmres')

    call lsolve_qmr(a, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=100)
    call check(out, xtrue, 'qmr')

    print '(a)', 'test_krylov: PASS'

contains

    subroutine check(out, expected, name)
        type(solver_result), intent(in) :: out !! Solver result to verify for convergence and numerical accuracy.
        real(dp), intent(in) :: expected(:) !! Expected solution vector.
        character(len=*), intent(in) :: name !! Short method name included in any failing error message.

        if (.not. out%converged) then
            print *, trim(name), ' flag=', out%flag, ' iter=', out%iter, ' err=', out%errors(size(out%errors))
            error stop 'krylov convergence'
        end if
        if (maxval(abs(out%x - expected)) > 2.0e-6_dp) then
            print *, trim(name), out%x
            error stop 'krylov solution'
        end if
    end subroutine check

end program test_krylov
