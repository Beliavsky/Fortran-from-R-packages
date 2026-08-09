program test_semidefinite
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: dmat(2,2), dvec(2), amat(2,2), bvec(2)
    type(dykstra_result) :: result

    dmat = reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2,2])
    dvec = [1.0_dp, 1.0_dp]
    amat = 0.0_dp
    amat(1,1) = 1.0_dp
    amat(2,2) = 1.0_dp
    bvec = 0.0_dp

    call dykstra_solve(dmat, dvec, amat, result, bvec=bvec)
    if (.not. result%converged) error stop "semidefinite problem did not converge"
    if (maxval(abs(result%solution - 0.5_dp)) > 2.0e-12_dp) then
        error stop "semidefinite regularized solution mismatch"
    end if
    if (abs(result%value + 0.5_dp) > 2.0e-12_dp) error stop "semidefinite objective mismatch"
end program test_semidefinite
