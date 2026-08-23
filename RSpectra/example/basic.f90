! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

program basic
    use rspectra
    implicit none
    real(dp) :: a(6,6), x(8,5)
    type(eigs_sym_result) :: er
    type(svds_result) :: sr
    integer :: i

    a = 0.0_dp
    do i = 1, 6
        a(i,i) = real(i,dp)
    end do
    er = eigs_sym(a, 2, 'LM')
    print '(a,2f12.6)', 'largest eigenvalues: ', er%values

    x = 0.0_dp
    x(1,1) = 5.0_dp
    x(2,2) = 3.0_dp
    x(3,3) = 1.0_dp
    sr = svds(x, 2)
    print '(a,2f12.6)', 'largest singular values: ', sr%d
end program basic
