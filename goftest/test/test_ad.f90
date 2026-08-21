program test_ad
    use goftest, only : dp, p_ad, q_ad, ad_statistic, ad_test_uniform, ad_inf_exact
    implicit none

    real(dp), parameter :: tol = 3.0e-11_dp
    real(dp) :: u(6), stat, pval, q

    call check_close(p_ad(0.2_dp), 0.009586709322645643_dp, tol, 'fast asymptotic cdf')
    call check_close(p_ad(0.5_dp), 0.25318227234000723_dp, tol, 'fast asymptotic cdf 2')
    call check_close(p_ad(1.0_dp), 0.6427139632681818_dp, tol, 'fast asymptotic cdf 3')
    call check_close(ad_inf_exact(0.5_dp), 0.25318562646965503_dp, 3.0e-12_dp, 'exact asymptotic cdf')
    call check_close(ad_inf_exact(2.0_dp), 0.9081632250587462_dp, 3.0e-12_dp, 'exact asymptotic cdf 2')
    call check_close(p_ad(0.5_dp, n=10), 0.2573659942339062_dp, tol, 'finite n cdf')
    call check_close(p_ad(2.0_dp, n=10), 0.9069353492122374_dp, tol, 'finite n cdf 2')

    q = q_ad(0.5_dp)
    call check_close(q, 0.7742327955823599_dp, 2.0e-10_dp, 'fast quantile')
    call check_close(p_ad(q), 0.5_dp, 2.0e-12_dp, 'fast inversion')
    q = q_ad(0.9_dp, fast=.false.)
    call check_close(q, 1.9329578327416093_dp, 2.0e-10_dp, 'exact quantile')
    call check_close(p_ad(q, fast=.false.), 0.9_dp, 2.0e-12_dp, 'exact inversion')
    q = q_ad(0.5_dp, n=20)
    call check_close(q, 0.7718581377639052_dp, 3.0e-10_dp, 'finite quantile')

    u = [0.07_dp, 0.18_dp, 0.32_dp, 0.58_dp, 0.76_dp, 0.91_dp]
    call check_close(ad_statistic(u), 0.18867283336985619_dp, 2.0e-14_dp, 'AD statistic')
    call ad_test_uniform(u, stat, pval)
    call check_close(stat, 0.18867283336985619_dp, 2.0e-14_dp, 'AD test statistic')
    call check_close(pval, 1.0_dp - p_ad(stat, n=6), 2.0e-14_dp, 'AD p-value')

    ! The Marsaglia finite-n correction can become slightly negative in the far lower tail.
    ! The Fortran public CDF clamps that invalid probability to zero.
    call check_close(p_ad(0.1_dp, n=10), 0.0_dp, 0.0_dp, 'finite-n lower-tail clamp')

    print '(a)', 'test_ad: PASS'

contains

    subroutine check_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual, expected, atol
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > atol) then
            write(*,'(a,2(1x,es24.16))') trim(label)//' mismatch:', actual, expected
            error stop 1
        end if
    end subroutine check_close

end program test_ad
