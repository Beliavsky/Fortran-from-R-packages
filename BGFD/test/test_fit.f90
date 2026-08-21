program test_fit
    use bgfd
    implicit none
    integer, parameter :: n = 500
    real(dp) :: x(n), p
    integer :: i, failures
    type(bgfd_fit_result) :: fit

    failures = 0
    do i = 1, n
        p = (real(i,dp)-0.5_dp)/real(n,dp)
        x(i) = q_bell_e(p, 1.2_dp, 0.7_dp)
    end do
    call m_bell_e(x, 1.0_dp, 0.9_dp, fit, method='B', max_iter=3000)
    if (fit%convergence /= 0) then
        print '(a,i0)', 'Bell-E convergence FAIL: ', fit%convergence
        failures = failures + 1
    end if
    call check_close('Bell-E alpha', fit%params(1), 1.2_dp, 0.04_dp, failures)
    call check_close('Bell-E lambda', fit%params(2), 0.7_dp, 0.09_dp, failures)
    if (.not. fit%pdf_integral_ok .or. .not. fit%cdf_endpoints_ok) failures = failures + 1

    do i = 1, n
        p = (real(i,dp)-0.5_dp)/real(n,dp)
        x(i) = q_cbell_w(p, 0.9_dp, 1.6_dp, 0.7_dp)
    end do
    call m_cbell_w(x, 0.8_dp, 1.4_dp, 0.9_dp, fit, method='B', max_iter=3000)
    if (fit%convergence /= 0) then
        print '(a,i0)', 'CBell-W convergence FAIL: ', fit%convergence
        failures = failures + 1
    end if
    call check_close('CBell-W alpha', fit%params(1), 0.9_dp, 0.06_dp, failures)
    call check_close('CBell-W beta', fit%params(2), 1.6_dp, 0.06_dp, failures)
    call check_close('CBell-W lambda', fit%params(3), 0.7_dp, 0.14_dp, failures)

    if (failures /= 0) error stop 1
    print '(a)', 'test_fit: PASS'
contains
    subroutine check_close(name, got, expected, atol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, expected, atol
        integer, intent(inout) :: failures
        if (abs(got-expected) > atol) then
            print '(a,2f14.7)', trim(name)//' FAIL: ', got, expected
            failures = failures + 1
        end if
    end subroutine check_close
end program test_fit
