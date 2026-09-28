program test_aux_stationary
    use rlinsolve, only : dp, sparse_matrix, solver_result, aux_fisch, &
        lsolve_jacobi, lsolve_gs, lsolve_sor, lsolve_ssor, sparse_to_dense
    implicit none

    real(dp), allocatable :: a(:, :)
    real(dp), allocatable :: as_dense(:, :)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: xtrue(:)
    type(sparse_matrix) :: as
    type(solver_result) :: out
    integer :: n

    call aux_fisch(3, a)
    call aux_fisch(3, as)
    if (size(a, 1) /= 9 .or. size(a, 2) /= 9) error stop "aux dense shape"
    if (as%nnz() /= 33) error stop "aux sparse nnz"
    as_dense = sparse_to_dense(as)
    if (maxval(abs(a - as_dense)) > 0.0_dp) error stop "aux sparse parity"

    n = size(a, 1)
    allocate(xtrue(n), b(n))
    xtrue = 1.0_dp
    b = matmul(a, xtrue)

    call lsolve_jacobi(a, b, out, xinit=0.0_dp*xtrue, reltol=1.0e-9_dp, maxiter=5000)
    if (.not. out%converged) error stop "jacobi convergence"
    if (maxval(abs(out%x - xtrue)) > 1.0e-6_dp) error stop "jacobi solution"

    call lsolve_gs(a, b, out, xinit=0.0_dp*xtrue, reltol=1.0e-10_dp, maxiter=2000)
    if (.not. out%converged) error stop "gs convergence"
    if (maxval(abs(out%x - xtrue)) > 1.0e-7_dp) error stop "gs solution"

    call lsolve_sor(as, b, out, xinit=0.0_dp*xtrue, reltol=1.0e-10_dp, maxiter=2000, weight=1.1_dp)
    if (.not. out%converged) error stop "sor convergence"
    if (maxval(abs(out%x - xtrue)) > 1.0e-7_dp) error stop "sor solution"

    call lsolve_ssor(as, b, out, xinit=0.0_dp*xtrue, reltol=1.0e-9_dp, maxiter=2000, weight=1.0_dp)
    if (.not. out%converged) error stop "ssor convergence"
    if (maxval(abs(out%x - xtrue)) > 1.0e-6_dp) error stop "ssor solution"

    print '(a)', 'test_aux_stationary: PASS'
end program test_aux_stationary
