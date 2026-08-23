program test_options
    use ecpdist, only: dp, pecp, qecp, secp, ecp_cumhaz
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    real(dp), parameter :: tol = 2.0e-12_dp
    real(dp) :: x, p
    integer :: fails

    fails = 0
    x = 0.73_dp
    p = 0.2_dp

    call check_close(pecp(x, 1.3_dp, 0.9_dp, 1.7_dp, lower_tail=.false.), &
        1.0_dp - pecp(x, 1.3_dp, 0.9_dp, 1.7_dp), tol, 'upper cdf', fails)
    call check_close(exp(pecp(x, 1.3_dp, 0.9_dp, 1.7_dp, log_p=.true.)), &
        pecp(x, 1.3_dp, 0.9_dp, 1.7_dp), tol, 'log cdf', fails)
    call check_close(qecp(p, 1.3_dp, 0.9_dp, 1.7_dp, lower_tail=.false.), &
        qecp(1.0_dp - p, 1.3_dp, 0.9_dp, 1.7_dp), tol, 'upper quantile', fails)
    call check_close(qecp(log(p), 1.3_dp, 0.9_dp, 1.7_dp, log_p=.true.), &
        qecp(p, 1.3_dp, 0.9_dp, 1.7_dp), tol, 'log probability quantile', fails)
    call check_close(secp(x, 1.3_dp, 0.9_dp, 1.7_dp, cum_haz=.true.), &
        ecp_cumhaz(x, 1.3_dp, 0.9_dp, 1.7_dp), tol, 'cumhaz', fails)
    call check_close(pecp(qecp(0.831_dp, 1.2_dp, 0.8_dp, -0.7_dp), &
        1.2_dp, 0.8_dp, -0.7_dp), 0.831_dp, 5.0e-12_dp, 'cdf quantile', fails)

    if (.not. ieee_is_nan(qecp(1.1_dp, 1.0_dp, 1.0_dp, 1.0_dp))) then
        print '(a)', 'invalid probability did not return NaN'
        fails = fails + 1
    end if

    if (fails /= 0) error stop 1
    print '(a)', 'test_options: PASS'

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

end program test_options
