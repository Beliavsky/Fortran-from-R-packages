! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program test_gen
    use rspectra
    implicit none
    real(dp) :: a(10,10)
    type(eigs_result) :: r
    type(eigs_opts) :: o
    integer :: fails
    real(dp) :: resid

    fails = 0
    a = 0.0_dp
    a(1,2) = -2.0_dp
    a(2,1) = 2.0_dp
    a(3,4) = -1.0_dp
    a(4,3) = 1.0_dp
    a(5,5) = 5.0_dp
    a(6,6) = 3.0_dp
    a(7,7) = -4.0_dp
    a(8,8) = -2.0_dp
    a(9,9) = 0.5_dp
    a(10,10) = 0.25_dp
    allocate(o%initvec(10))
    o%initvec = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp,9.0_dp,10.0_dp]

    r = eigs(a, 2, 'LR', opts=o)
    if (r%info /= 0 .or. r%nconv /= 2) fails = fails + 1
    if (abs(real(r%values(1),dp) - 5.0_dp) > 1.0e-8_dp) fails = fails + 1
    if (abs(real(r%values(2),dp) - 3.0_dp) > 1.0e-8_dp) fails = fails + 1
    resid = maxval(abs(matmul(cmplx(a,0.0_dp,dp),r%vectors) - &
        r%vectors * spread(r%values,1,10)))
    if (resid > 1.0e-7_dp) fails = fails + 1

    r = eigs(a, 2, 'LI', opts=o)
    if (abs(abs(aimag(r%values(1))) - 2.0_dp) > 1.0e-8_dp) fails = fails + 1
    if (abs(abs(aimag(r%values(2))) - 2.0_dp) > 1.0e-8_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_gen: PASS'
    else
        print '(a,i0)', 'test_gen: FAIL ', fails
        error stop 1
    end if
end program test_gen
