! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program test_sparse
    use rspectra
    implicit none
    type(csr_operator) :: op
    type(eigs_sym_result) :: r
    integer :: row_ptr(11), col_ind(10), i, fails
    real(dp) :: val(10)

    fails = 0
    do i = 1, 11
        row_ptr(i) = i
    end do
    do i = 1, 10
        col_ind(i) = i
        val(i) = real(i,dp)
    end do
    op = make_csr_operator(10,10,row_ptr,col_ind,val)
    r = eigs_sym(op, 3, 'LM')
    if (r%info /= 0) fails = fails + 1
    if (maxval(abs(r%values - [10.0_dp,9.0_dp,8.0_dp])) > 1.0e-8_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_sparse: PASS'
    else
        print '(a,i0)', 'test_sparse: FAIL ', fails
        error stop 1
    end if
end program test_sparse
