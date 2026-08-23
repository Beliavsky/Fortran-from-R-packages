! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program test_shift
    use rspectra
    implicit none
    real(dp) :: a(8,8), g(8,8)
    type(eigs_sym_result) :: rs
    type(eigs_result) :: rg
    integer :: i, fails

    fails = 0
    a = 0.0_dp
    a(1,1) = -4.0_dp
    a(2,2) = -1.0_dp
    a(3,3) = 0.2_dp
    a(4,4) = 2.0_dp
    a(5,5) = 5.0_dp
    a(6,6) = 8.0_dp
    a(7,7) = 10.0_dp
    a(8,8) = 15.0_dp
    rs = eigs_sym(a, 2, 'LM', sigma=0.0_dp)
    if (maxval(abs(rs%values - [0.2_dp,-1.0_dp])) > 1.0e-7_dp) fails = fails + 1

    g = 0.0_dp
    g(1,2) = -2.0_dp
    g(2,1) = 2.0_dp
    g(3,4) = -1.0_dp
    g(4,3) = 1.0_dp
    do i = 5, 8
        g(i,i) = real(i-4,dp)
    end do
    rg = eigs(g, 2, 'LM', sigma=cmplx(0.0_dp,1.8_dp,dp))
    if (abs(abs(aimag(rg%values(1))) - 2.0_dp) > 1.0e-10_dp) fails = fails + 1
    if (abs(abs(aimag(rg%values(2))) - 1.0_dp) > 1.0e-10_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_shift: PASS'
    else
        print '(a,i0)', 'test_shift: FAIL ', fails
        error stop 1
    end if
end program test_shift
