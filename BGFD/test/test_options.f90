program test_options
    use bgfd
    implicit none
    integer :: failures
    real(dp) :: x, p, q, d

    failures = 0
    x = 0.8_dp
    p = p_bell_e(x, 1.2_dp, 0.7_dp)
    call check_close('upper tail', p_bell_e(x,1.2_dp,0.7_dp,lower_tail=.false.), &
        1.0_dp-p, 1.0e-13_dp, failures)
    call check_close('log p', exp(p_bell_e(x,1.2_dp,0.7_dp,log_p=.true.)), &
        p, 1.0e-13_dp, failures)
    q = q_bell_e(log(0.3_dp),1.2_dp,0.7_dp,log_p=.true.)
    call check_close('log quantile p', p_bell_e(q,1.2_dp,0.7_dp), 0.3_dp, 1.0e-12_dp, failures)
    q = q_bell_e(0.7_dp,1.2_dp,0.7_dp,lower_tail=.false.)
    call check_close('upper quantile', p_bell_e(q,1.2_dp,0.7_dp), 0.3_dp, 1.0e-12_dp, failures)
    d = d_cbell_b(x,1.1_dp,2.2_dp,1.4_dp,0.7_dp)
    call check_close('log density', exp(d_cbell_b(x,1.1_dp,2.2_dp,1.4_dp,0.7_dp,.true.)), &
        d, 1.0e-12_dp, failures)
    call check_close('survival', s_cbell_b(x,1.1_dp,2.2_dp,1.4_dp,0.7_dp), &
        1.0_dp-p_cbell_b(x,1.1_dp,2.2_dp,1.4_dp,0.7_dp), 1.0e-13_dp, failures)

    if (failures /= 0) error stop 1
    print '(a)', 'test_options: PASS'
contains
    subroutine check_close(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, expected, tol
        integer, intent(inout) :: failures
        if (abs(got-expected) > tol*max(1.0_dp,abs(expected))) then
            print '(a,2es24.15)', trim(name)//' FAIL: ', got, expected
            failures = failures + 1
        end if
    end subroutine check_close
end program test_options
