program test_status
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: dmat(2,2), dvec(2), amat(2,2), bvec(2)
    type(dykstra_result) :: result

    dmat = 0.0_dp
    dmat(1,1) = 1.0_dp
    dmat(2,2) = 1.0_dp
    dvec = [-1.0_dp, -1.0_dp]
    amat = 0.0_dp
    amat(1,1) = 1.0_dp
    amat(2,2) = 1.0_dp
    bvec = 0.0_dp

    call dykstra_solve(dmat, dvec, amat, result, bvec=bvec, maxit=1, eps=0.0_dp)
    if (result%converged) error stop "maxit=1 should stop before convergence"
    if (result%status /= 1) error stop "unexpected nonconvergence status"
    if (result%iterations /= 1) error stop "iteration count mismatch"
end program test_status
