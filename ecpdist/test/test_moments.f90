program test_moments
    use ecpdist, only: dp, ecp_integral_result, ecp_kmoment, ecp_kmoment_cond, ecp_mrl
    implicit none
    type(ecp_integral_result) :: m1, m2, c1, r
    integer :: fails

    fails = 0
    m1 = ecp_kmoment(1, 0.1_dp, 0.5_dp, -0.2_dp)
    m2 = ecp_kmoment(2, 0.1_dp, 0.5_dp, -0.2_dp)
    c1 = ecp_kmoment_cond(5.0_dp, 1, 0.1_dp, 0.5_dp, -0.2_dp)
    r = ecp_mrl(5.0_dp, 0.1_dp, 0.5_dp, -0.2_dp)

    call check_status(m1%status, 'moment 1 status', fails)
    call check_status(m2%status, 'moment 2 status', fails)
    call check_status(c1%status, 'conditional status', fails)
    call check_status(r%status, 'mrl status', fails)

    call check_close(m1%estimate, 4.6852446008722524_dp, 2.0e-8_dp, 'moment 1', fails)
    call check_close(m2%estimate, 35.135290802225067_dp, 2.0e-7_dp, 'moment 2', fails)
    call check_close(c1%estimate, 8.3043344145801956_dp, 3.0e-8_dp, 'conditional mean', fails)
    call check_close(r%estimate, 3.3043344145801956_dp, 3.0e-8_dp, 'mrl', fails)

    c1 = ecp_kmoment_cond(0.0_dp, 1, 0.1_dp, 0.5_dp, -0.2_dp)
    call check_close(c1%estimate, m1%estimate, 3.0e-9_dp, 'conditional at zero', fails)

    if (fails /= 0) error stop 1
    print '(a)', 'test_moments: PASS'

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

    subroutine check_status(status, label, nfail)
        integer, intent(in) :: status
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (status /= 0) then
            print '(a,i0)', trim(label)//': ', status
            nfail = nfail + 1
        end if
    end subroutine check_status

end program test_moments
