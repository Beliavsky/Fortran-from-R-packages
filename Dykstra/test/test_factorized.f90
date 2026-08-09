program test_factorized
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: rinv(2,2), dvec(2), amat(2,3), bvec(3)
    type(dykstra_result) :: result

    rinv = 0.0_dp
    rinv(1,1) = 0.5_dp
    rinv(2,2) = 1.0_dp / 3.0_dp
    dvec = [4.0_dp, 9.0_dp]
    amat = 0.0_dp
    amat(1,1) = 1.0_dp
    amat(2,2) = 1.0_dp
    amat(:,3) = 1.0_dp
    bvec = [0.0_dp, 0.0_dp, 1.5_dp]

    call dykstra_solve(rinv, dvec, amat, result, bvec=bvec, factorized=.true.)
    if (.not. result%converged) error stop "factorized problem did not converge"
    if (maxval(abs(result%solution - [1.0_dp, 1.0_dp])) > 5.0e-13_dp) then
        error stop "factorized solution mismatch"
    end if
    if (.not. ieee_is_nan(result%value)) error stop "factorized objective must be NaN"

    ! Dykstra 1.0-0 reports this quantity using dvec / Rinv**2 on diagonal inputs.
    if (maxval(abs(result%unconstrained - [16.0_dp, 81.0_dp])) > 5.0e-13_dp) then
        error stop "factorized source-compatible unconstrained output mismatch"
    end if
end program test_factorized
