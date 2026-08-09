program simplex_qp
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: dmat(3,3), dvec(3), amat(3,4), bvec(4)
    type(dykstra_result) :: result

    dmat = 0.0_dp
    dmat(1,1) = 1.0_dp
    dmat(2,2) = 1.0_dp
    dmat(3,3) = 1.0_dp
    dvec = [1.0_dp, 1.5_dp, 1.0_dp]
    amat = 0.0_dp
    amat(:,1) = 1.0_dp
    amat(1,2) = 1.0_dp
    amat(2,3) = 1.0_dp
    amat(3,4) = 1.0_dp
    bvec = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

    call dykstra_solve(dmat, dvec, amat, result, bvec=bvec, meq=1)
    print '(a,3f14.8)', 'solution: ', result%solution
    print '(a,es16.8)', 'objective: ', result%value
    print '(a,i0)', 'cycles: ', result%iterations
    print '(a,l1)', 'converged: ', result%converged
end program simplex_qp
