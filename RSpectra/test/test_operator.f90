! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module test_operator_mod
    use rspectra
    implicit none
    type, extends(linear_operator) :: diag_operator
        real(dp), allocatable :: d(:)
    contains
        procedure :: prod => diag_prod
        procedure :: tprod => diag_tprod
    end type diag_operator
contains
    subroutine diag_prod(self, x, y)
        class(diag_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = self%d * x
    end subroutine diag_prod
    subroutine diag_tprod(self, x, y)
        class(diag_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = self%d * x
    end subroutine diag_tprod
end module test_operator_mod

program test_operator
    use test_operator_mod
    implicit none
    type(diag_operator) :: op
    type(eigs_sym_result) :: er
    type(svds_result) :: sr
    integer :: i, fails

    fails = 0
    op%nrow = 12
    op%ncol = 12
    allocate(op%d(12))
    op%d = [(real(i,dp), i=1,12)]
    er = eigs_sym(op, 3, 'LM')
    if (maxval(abs(er%values - [12.0_dp,11.0_dp,10.0_dp])) > 1.0e-8_dp) fails = fails + 1
    sr = svds(op, 3)
    if (maxval(abs(sr%d - [12.0_dp,11.0_dp,10.0_dp])) > 1.0e-8_dp) fails = fails + 1

    if (fails == 0) then
        print '(a)', 'test_operator: PASS'
    else
        print '(a,i0)', 'test_operator: FAIL ', fails
        error stop 1
    end if
end program test_operator
