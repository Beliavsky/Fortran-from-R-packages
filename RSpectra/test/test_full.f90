! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program test_full
    use rspectra
    implicit none
    real(dp) :: a(4,4), s(4,4)
    type(eigs_result) :: rg
    type(eigs_sym_result) :: rs
    type(svds_result) :: rv
    integer :: i, fails
    fails = 0
    a = 0.0_dp
    a(1,2) = -1.0_dp
    a(2,1) = 1.0_dp
    a(3,3) = 3.0_dp
    a(4,4) = -2.0_dp
    rg = eigs(a,4,'LM')
    if (rg%nconv /= 4 .or. rg%info /= 0) fails = fails + 1
    if (abs(abs(rg%values(1)) - 3.0_dp) > 1.0e-10_dp) fails = fails + 1

    s = 0.0_dp
    do i = 1, 4
        s(i,i) = real(i,dp)
    end do
    rs = eigs_sym(s,4,'LA')
    if (maxval(abs(rs%values - [4.0_dp,3.0_dp,2.0_dp,1.0_dp])) > 1.0e-10_dp) fails = fails + 1
    rv = svds(s,4)
    if (maxval(abs(rv%d - [4.0_dp,3.0_dp,2.0_dp,1.0_dp])) > 1.0e-10_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_full: PASS'
    else
        print '(a,i0)', 'test_full: FAIL ', fails
        error stop 1
    end if
end program test_full
