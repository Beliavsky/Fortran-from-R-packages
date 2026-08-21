program test_cvm
    use goftest, only : dp, p_cvm, q_cvm, cvm_statistic, cvm_test_uniform
    implicit none

    real(dp), parameter :: tol = 2.0e-10_dp
    real(dp) :: u(6), stat, pval, q

    call check_close(p_cvm(0.03_dp), 0.02383154944170798_dp, tol, 'asymptotic cdf')
    call check_close(p_cvm(0.10_dp), 0.41512656159320305_dp, tol, 'asymptotic cdf 2')
    call check_close(p_cvm(0.50_dp), 0.9601667824343911_dp, tol, 'asymptotic cdf 3')
    call check_close(p_cvm(0.10_dp, n=10), 0.4060001239613563_dp, tol, 'finite cdf')
    call check_close(p_cvm(0.50_dp, n=100), 0.9603805253156238_dp, tol, 'finite cdf 2')

    q = q_cvm(0.5_dp)
    call check_close(q, 0.11887955098036691_dp, 5.0e-10_dp, 'asymptotic quantile')
    call check_close(p_cvm(q), 0.5_dp, 5.0e-10_dp, 'asymptotic inversion')
    q = q_cvm(0.9_dp, n=10)
    call check_close(q, 0.34514342066466397_dp, 5.0e-10_dp, 'finite quantile')
    call check_close(p_cvm(q, n=10), 0.9_dp, 5.0e-10_dp, 'finite inversion')

    u = [0.07_dp, 0.18_dp, 0.32_dp, 0.58_dp, 0.76_dp, 0.91_dp]
    call check_close(cvm_statistic(u), 0.028466666666666668_dp, 2.0e-15_dp, 'CvM statistic')
    call cvm_test_uniform(u, stat, pval)
    call check_close(stat, 0.028466666666666668_dp, 2.0e-15_dp, 'CvM test statistic')
    call check_close(pval, p_cvm(stat, n=6, lower_tail=.false.), 2.0e-14_dp, 'CvM p-value')

    print '(a)', 'test_cvm: PASS'

contains

    subroutine check_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual, expected, atol
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > atol) then
            write(*,'(a,2(1x,es24.16))') trim(label)//' mismatch:', actual, expected
            error stop 1
        end if
    end subroutine check_close

end program test_cvm
