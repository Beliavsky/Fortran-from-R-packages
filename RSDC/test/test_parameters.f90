! SPDX-License-Identifier: GPL-3.0-only
program test_parameters
    use rsdc, only: dp, partial_to_correlation, correlation_to_partial
    use rsdc, only: transition_at, fixed_transition_from_parameters
    implicit none
    real(dp) :: z(6), z2(6), r(4,4), p(3,3), p2(2,2), beta(3,4), x(2)
    logical :: ok
    z = [0.2_dp, -0.1_dp, 0.3_dp, 0.15_dp, -0.25_dp, 0.4_dp]
    call partial_to_correlation(z, 4, r)
    call correlation_to_partial(r, z2)
    call check(maxval(abs(z-z2)) < 1.0e-10_dp, 'partial round trip')
    call check(maxval(abs(r-transpose(r))) < 1.0e-12_dp, 'correlation symmetry')
    beta = reshape([0.2_dp,0.1_dp,-0.2_dp,0.3_dp, 0.4_dp,-0.1_dp,0.2_dp,0.1_dp, &
                    -0.3_dp,0.2_dp,0.1_dp,-0.2_dp], shape(beta), order=[2,1])
    x = [1.0_dp, 0.5_dp]
    call transition_at(beta, x, p)
    call check(maxval(abs(sum(p,dim=2)-1.0_dp)) < 1.0e-12_dp, 'softmax rows')
    call check(all(p > 0.0_dp), 'softmax positivity')
    call fixed_transition_from_parameters([0.8_dp,0.7_dp], 2, p2, ok)
    call check(ok .and. abs(p2(1,2)-0.2_dp)<1.0e-12_dp, 'fixed transition')
    print '(a)', 'test_parameters: PASS'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) error stop message
    end subroutine check
end program test_parameters
