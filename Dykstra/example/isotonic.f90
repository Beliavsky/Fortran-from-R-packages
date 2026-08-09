program isotonic
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: dmat(5,5), y(5), amat(5,4)
    type(dykstra_result) :: result
    integer :: i

    dmat = 0.0_dp
    do i = 1, 5
        dmat(i,i) = 1.0_dp
    end do
    y = [3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 4.0_dp]
    amat = 0.0_dp
    do i = 1, 4
        amat(i,i) = -1.0_dp
        amat(i+1,i) = 1.0_dp
    end do

    call dykstra_solve(dmat, y, amat, result)
    print '(a,5f10.4)', 'data:      ', y
    print '(a,5f10.4)', 'isotonic:  ', result%solution
end program isotonic
