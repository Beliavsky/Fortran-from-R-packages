program test_isotonic
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: dmat(3,3), y(3), amat(3,2), bvec(2)
    type(dykstra_result) :: result

    dmat = 0.0_dp
    dmat(1,1) = 1.0_dp
    dmat(2,2) = 1.0_dp
    dmat(3,3) = 1.0_dp
    y = [3.0_dp, 1.0_dp, 2.0_dp]
    amat(:,1) = [-1.0_dp, 1.0_dp, 0.0_dp]
    amat(:,2) = [0.0_dp, -1.0_dp, 1.0_dp]
    bvec = 0.0_dp

    call dykstra_solve(dmat, y, amat, result, bvec=bvec)
    if (.not. result%converged) error stop "isotonic regression did not converge"
    if (maxval(abs(result%solution - 2.0_dp)) > 5.0e-13_dp) then
        error stop "isotonic solution mismatch"
    end if
end program test_isotonic
