program test_pmultinom
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
    use pmultinom_module, only : dp, pmultinom
    implicit none

    integer :: failures
    real(dp) :: inf, ninf

    failures = 0
    inf = ieee_value(1.0_dp, ieee_positive_inf)
    ninf = ieee_value(1.0_dp, ieee_negative_inf)

    call check_close('all upper zero, n=0', pmultinom(0, [0.5_dp, 0.5_dp], upper=[0.0_dp, 0.0_dp]), 1.0_dp, 1.0e-13_dp)
    call check_close('all upper zero, n=1', pmultinom(1, [0.5_dp, 0.5_dp], upper=[0.0_dp, 0.0_dp]), 0.0_dp, 1.0e-13_dp)
    call check_close('certain upper', pmultinom(10, [0.5_dp, 0.5_dp], upper=[10.0_dp, 10.0_dp]), 1.0_dp, 1.0e-12_dp)
    call check_close('certain lower', pmultinom(10, [0.5_dp, 0.5_dp], lower=[-1.0_dp, -1.0_dp]), 1.0_dp, 1.0e-12_dp)
    call check_close('unbounded', pmultinom(10, [0.5_dp, 0.5_dp], lower=[ninf, ninf], upper=[inf, inf]), 1.0_dp, 1.0e-12_dp)
    call check_close('one category upper fail', pmultinom(10, [1.0_dp], upper=[5.0_dp]), 0.0_dp, 0.0_dp)
    call check_close('one category upper pass', pmultinom(5, [1.0_dp], upper=[5.0_dp]), 1.0_dp, 0.0_dp)
    call check_close('one category lower pass', pmultinom(10, [1.0_dp], lower=[5.0_dp]), 1.0_dp, 0.0_dp)
    call check_close('zero-prob lower fails', pmultinom(5, [1.0_dp, 0.0_dp], lower=[0.0_dp, 0.0_dp]), 0.0_dp, 1.0e-13_dp)
    call check_close('zero-prob allowed', pmultinom(5, [1.0_dp, 0.0_dp], lower=[0.0_dp, -1.0_dp]), 1.0_dp, 1.0e-12_dp)

    ! P(N1 <= 3, N2 <= 5) for n=5 and p=(1/2,1/2) = P(N1 >= 0 and N1 <= 3) = 26/32.
    call check_close('small exact cdf', pmultinom(5, [0.5_dp, 0.5_dp], upper=[3.0_dp, 5.0_dp]), 26.0_dp/32.0_dp, 5.0e-12_dp)

    ! P(N1 > 2, N2 > -1) = P(N1 >= 3) for Binomial(5, 1/2) = 16/32.
    call check_close('small exact tail', pmultinom(5, [0.5_dp, 0.5_dp], lower=[2.0_dp, -1.0_dp]), 0.5_dp, 5.0e-12_dp)

    ! General two-sided event: 1 < N1 <= 3 and -1 < N2 <= 5 means N1 = 2 or 3.
    call check_close('small between', pmultinom(5, [0.5_dp, 0.5_dp], lower=[1.0_dp, -1.0_dp], &
                     upper=[3.0_dp, 5.0_dp]), 20.0_dp/32.0_dp, 5.0e-12_dp)

    ! Upstream large-size test: equivalent to dbinom(10, 100, 0.5).
    call check_close('large size binomial reduction', pmultinom(100, [0.5_dp, 0.5_dp], upper=[10.0_dp, 0.0_dp]), &
                     binom_pmf(10, 100, 0.5_dp), 1.0e-11_dp)

    call test_large_category_count()
    call check_close('Wolf fair-die cdf', pmultinom(20000, spread(1.0_dp/6.0_dp, 1, 6), &
                     upper=spread(3630.0_dp, 1, 6)), 1.0_dp - 7.379909e-8_dp, 5.0e-11_dp)

    if (failures /= 0) then
        print '(a,i0)', 'test_pmultinom: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_pmultinom: PASS'

contains

    subroutine check_close(label, got, expected, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: got, expected, tol
        real(dp) :: scale

        scale = max(1.0_dp, abs(expected))
        if (abs(got - expected) > tol * scale) then
            failures = failures + 1
            print '(a,2(1x,es24.16))', trim(label)//' got/expected:', got, expected
        end if
    end subroutine check_close


    subroutine test_large_category_count()
        real(dp), allocatable :: p(:), u(:)
        real(dp) :: got

        allocate(p(10000), u(10000))
        p = 1.0e-4_dp
        u = 0.0_dp
        u(1) = 1.0_dp
        got = pmultinom(1, p, upper=u)
        call check_close('10000 categories', got, 1.0e-4_dp, 1.0e-11_dp)
    end subroutine test_large_category_count

    real(dp) function binom_pmf(k, n, p) result(ans)
        integer, intent(in) :: k, n
        real(dp), intent(in) :: p
        ans = exp(log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) - &
                  log_gamma(real(n - k + 1, dp)) + real(k, dp) * log(p) + &
                  real(n-k, dp) * log(1.0_dp-p))
    end function binom_pmf

end program test_pmultinom
