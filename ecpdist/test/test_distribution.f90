program test_distribution
    use ecpdist, only: dp, decp, pecp, qecp, secp, hecp
    implicit none
    real(dp), parameter :: tol = 5.0e-12_dp
    real(dp) :: h, deriv, x
    integer :: fails

    fails = 0
    call check_close(decp(1.0_dp, 1.2_dp, 0.8_dp, -0.7_dp), &
        0.25056194470311317_dp, tol, 'negative-phi density', fails)
    call check_close(pecp(1.0_dp, 1.2_dp, 0.8_dp, -0.7_dp), &
        0.90813353915131078_dp, tol, 'negative-phi cdf', fails)
    call check_close(qecp(0.5_dp, 1.2_dp, 0.8_dp, -0.7_dp), &
        0.28719370802889544_dp, tol, 'negative-phi quantile', fails)

    call check_close(decp(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp), &
        0.64469318224302199_dp, tol, 'density reference', fails)
    call check_close(pecp(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp), &
        0.74022864902429551_dp, tol, 'cdf reference', fails)
    call check_close(qecp(0.9_dp, 1.0_dp, 1.0_dp, 1.0_dp), &
        1.3160834424264862_dp, tol, 'quantile reference', fails)

    x = 0.8_dp
    h = 1.0e-5_dp
    deriv = (pecp(x + h, 1.1_dp, 0.7_dp, 2.0_dp) - &
        pecp(x - h, 1.1_dp, 0.7_dp, 2.0_dp))/(2.0_dp*h)
    call check_close(deriv, decp(x, 1.1_dp, 0.7_dp, 2.0_dp), 2.0e-9_dp, &
        'cdf derivative', fails)
    call check_close(hecp(x, 1.1_dp, 0.7_dp, 2.0_dp), &
        decp(x, 1.1_dp, 0.7_dp, 2.0_dp)/secp(x, 1.1_dp, 0.7_dp, 2.0_dp), &
        2.0e-12_dp, 'hazard identity', fails)

    if (fails /= 0) error stop 1
    print '(a)', 'test_distribution: PASS'

contains

    subroutine check_close(got, expected, atol, label, nfail)
        real(dp), intent(in) :: got, expected, atol
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (abs(got - expected) > atol) then
            print '(a,2es24.15)', trim(label)//': ', got, expected
            nfail = nfail + 1
        end if
    end subroutine check_close

end program test_distribution
