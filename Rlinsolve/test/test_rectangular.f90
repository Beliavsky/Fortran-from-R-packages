program test_rectangular
    use rlinsolve, only : dp, solver_result, lsolve_cg, lsolve_jacobi
    implicit none

    real(dp) :: a(5, 2)
    real(dp) :: b(5)
    real(dp) :: x0(2)
    real(dp) :: xtrue(2)
    type(solver_result) :: out

    a = reshape([1.0_dp, 0.0_dp, &
                 0.0_dp, 1.0_dp, &
                 1.0_dp, 1.0_dp, &
                 2.0_dp, 1.0_dp, &
                 1.0_dp, 2.0_dp], [5, 2])
    xtrue = [2.0_dp, -1.0_dp]
    b = matmul(a, xtrue)
    x0 = 0.0_dp

    call lsolve_cg(a, b, out, xinit=x0, reltol=1.0e-12_dp, maxiter=100)
    if (.not. out%converged) error stop "rectangular cg convergence"
    if (maxval(abs(out%x - xtrue)) > 1.0e-8_dp) error stop "rectangular cg solution"

    call lsolve_jacobi(a, b, out, xinit=x0, reltol=1.0e-9_dp, maxiter=5000)
    if (.not. out%converged) error stop "rectangular jacobi convergence"
    if (maxval(abs(out%x - xtrue)) > 1.0e-6_dp) error stop "rectangular jacobi solution"

    print '(a)', 'test_rectangular: PASS'
end program test_rectangular
