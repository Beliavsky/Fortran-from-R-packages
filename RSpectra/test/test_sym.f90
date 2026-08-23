! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program test_sym
    use rspectra
    implicit none
    real(dp) :: a(8,8)
    type(eigs_sym_result) :: r
    type(eigs_opts) :: o
    integer :: i, fails
    real(dp) :: err

    fails = 0
    a = 0.0_dp
    do i = 1, 8
        a(i,i) = real(i,dp)
    end do
    allocate(o%initvec(8))
    o%initvec = [(real(i,dp), i=1,8)]

    r = eigs_sym(a, 3, 'LM', opts=o)
    if (r%info /= 0 .or. r%nconv /= 3) fails = fails + 1
    if (maxval(abs(r%values - [8.0_dp,7.0_dp,6.0_dp])) > 1.0e-9_dp) fails = fails + 1
    err = maxval(abs(matmul(a,r%vectors) - r%vectors * spread(r%values,1,8)))
    if (err > 1.0e-8_dp) fails = fails + 1

    r = eigs_sym(a, 3, 'SM', opts=o)
    if (maxval(abs(r%values - [1.0_dp,2.0_dp,3.0_dp])) > 1.0e-9_dp) fails = fails + 1

    r = eigs_sym(a, 4, 'BE', opts=o)
    if (maxval(abs(r%values - [8.0_dp,1.0_dp,7.0_dp,2.0_dp])) > 1.0e-8_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_sym: PASS'
    else
        print '(a,i0)', 'test_sym: FAIL ', fails
        error stop 1
    end if
end program test_sym
