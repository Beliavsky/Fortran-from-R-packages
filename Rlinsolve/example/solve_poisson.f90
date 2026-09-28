program solve_poisson
    use rlinsolve, only : dp, solver_result, aux_fisch, lsolve_cg
    implicit none

    real(dp), allocatable :: a(:, :)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: x0(:)
    real(dp), allocatable :: xtrue(:)
    type(solver_result) :: out

    call aux_fisch(4, a)
    allocate(xtrue(size(a, 1)), b(size(a, 1)), x0(size(a, 1)))
    xtrue = 1.0_dp
    b = matmul(a, xtrue)
    x0 = 0.0_dp

    call lsolve_cg(a, b, out, xinit=x0, reltol=1.0e-10_dp, maxiter=200)
    print '(a,l1)', 'converged = ', out%converged
    print '(a,es12.4)', 'max abs error = ', maxval(abs(out%x - xtrue))
end program solve_poisson
