! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program test_svds
    use rspectra
    implicit none
    real(dp) :: a(8,5), b(8,5)
    type(svds_result) :: r, ref
    integer :: fails

    fails = 0
    a = 0.0_dp
    a(1,1) = 5.0_dp
    a(2,2) = 3.0_dp
    a(3,3) = 1.0_dp
    a(4,4) = 0.5_dp
    a(5,5) = 0.2_dp
    r = svds(a,3)
    if (r%info /= 0 .or. r%nconv /= 3) fails = fails + 1
    if (maxval(abs(r%d - [5.0_dp,3.0_dp,1.0_dp])) > 1.0e-8_dp) fails = fails + 1
    if (maxval(abs(matmul(a,r%v) - r%u * spread(r%d,1,8))) > 1.0e-7_dp) fails = fails + 1

    r = svds(transpose(a),3)
    if (r%info /= 0 .or. r%nconv /= 3) fails = fails + 1
    if (maxval(abs(r%d - [5.0_dp,3.0_dp,1.0_dp])) > 1.0e-8_dp) fails = fails + 1
    if (maxval(abs(matmul(transpose(a),r%v) - r%u * spread(r%d,1,5))) > 1.0e-7_dp) fails = fails + 1

    b = reshape([ &
        1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp, &
        2.0_dp,1.0_dp,4.0_dp,3.0_dp,6.0_dp,5.0_dp,8.0_dp,7.0_dp, &
        1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp, &
        3.0_dp,1.0_dp,2.0_dp,4.0_dp,2.0_dp,5.0_dp,1.0_dp,6.0_dp, &
        2.0_dp,5.0_dp,1.0_dp,3.0_dp,4.0_dp,2.0_dp,6.0_dp,1.0_dp], [8,5])
    r = svds(b,3,center=.true.,scale=.true.)
    ref = svds(b,5,center=.true.,scale=.true.)
    if (maxval(abs(r%d - ref%d(1:3))) > 1.0e-7_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_svds: PASS'
    else
        print '(a,i0)', 'test_svds: FAIL ', fails
        error stop 1
    end if
end program test_svds
