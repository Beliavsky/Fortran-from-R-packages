program test_rng
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use sadists
    implicit none
    integer, parameter :: n = 30000
    integer :: failures
    integer :: seed_size, i
    integer, allocatable :: seed(:)
    real(dp), allocatable :: x(:)
    real(dp) :: k(2), m(2), target, sd

    failures = 0
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    do i = 1, seed_size
        seed(i) = 104729 + 7919*i
    end do
    call random_seed(put=seed)
    allocate(x(n))

    call rdnbeta(x, 10.0_dp, 20.0_dp, 0.0_dp, 0.0_dp)
    target = 1.0_dp/3.0_dp
    sd = sqrt((5.0_dp*10.0_dp) / (15.0_dp**2 * 16.0_dp))
    call check_mean(x, target, sd, 'rdnbeta', failures)

    call rdneta(x, 20.0_dp, 0.0_dp, 0.0_dp)
    call check_mean(x, 0.0_dp, sqrt(1.0_dp/21.0_dp), 'rdneta', failures)

    call dnf_moments(40.0_dp, 80.0_dp, 1.5_dp, 0.0_dp, m)
    call rdnf(x, 40.0_dp, 80.0_dp, 1.5_dp, 0.0_dp)
    call check_mean(x, m(1), sqrt(max(0.0_dp, m(2)-m(1)**2)), &
        'rdnf', failures)

    call dnt_moments(75.0_dp, 2.0_dp, 0.0_dp, m)
    call rdnt(x, 75.0_dp, 2.0_dp, 0.0_dp)
    call check_mean(x, m(1), sqrt(max(0.0_dp, m(2)-m(1)**2)), &
        'rdnt', failures)

    call lambdap_cumulants(50.0_dp, 1.5_dp, k)
    call rlambdap(x, 50.0_dp, 1.5_dp)
    call check_mean(x, k(1), sqrt(k(2)), 'rlambdap', failures)

    call upsilon_cumulants([30.0_dp, 50.0_dp], [-0.5_dp, 1.0_dp], k)
    call rupsilon(x, [30.0_dp, 50.0_dp], [-0.5_dp, 1.0_dp])
    call check_mean(x, k(1), sqrt(k(2)), 'rupsilon', failures)

    call kprime_cumulants(50.0_dp, 80.0_dp, 0.5_dp, 1.5_dp, k)
    call rkprime(x, 50.0_dp, 80.0_dp, 0.5_dp, 1.5_dp)
    call check_mean(x, k(1), sqrt(k(2)), 'rkprime', failures)

    call sumchisqpow_cumulants([-1.0_dp, 1.0_dp], [100.0_dp, 200.0_dp], &
        [0.0_dp, 1.0_dp], [1.0_dp, 0.5_dp], k)
    call rsumchisqpow(x, [-1.0_dp, 1.0_dp], [100.0_dp, 200.0_dp], &
        [0.0_dp, 1.0_dp], [1.0_dp, 0.5_dp])
    call check_mean(x, k(1), sqrt(k(2)), 'rsumchisqpow', failures)

    call sumlogchisq_cumulants([1.0_dp, -0.5_dp], [20.0_dp, 30.0_dp], &
        [2.0_dp, 3.0_dp], k)
    call rsumlogchisq(x, [1.0_dp, -0.5_dp], [20.0_dp, 30.0_dp], &
        [2.0_dp, 3.0_dp])
    call check_mean(x, k(1), sqrt(k(2)), 'rsumlogchisq', failures)

    call rprodnormal(x, [1.2_dp, -0.7_dp], [0.5_dp, 0.3_dp])
    if (.not. all(ieee_is_finite(x))) then
        print *, 'rprodnormal generated non-finite values'
        failures = failures + 1
    end if
    call rprodchisqpow(x, [20.0_dp, 30.0_dp], [2.0_dp, 3.0_dp], &
        [1.0_dp, 0.5_dp])
    if (any(x <= 0.0_dp) .or. .not. all(ieee_is_finite(x))) then
        print *, 'rprodchisqpow support failure'
        failures = failures + 1
    end if
    call rproddnf(x, [10.0_dp, 20.0_dp], [100.0_dp, 80.0_dp], &
        [1.0_dp, 0.5_dp], [0.2_dp, 1.0_dp])
    if (any(x <= 0.0_dp) .or. .not. all(ieee_is_finite(x))) then
        print *, 'rproddnf support failure'
        failures = failures + 1
    end if

    if (failures == 0) then
        print *, 'test_rng: PASS'
    else
        print *, 'test_rng: FAIL', failures
        error stop 1
    end if

contains

    subroutine check_mean(y, expected, sigma, name, failures)
        real(dp), intent(in) :: y(:), expected, sigma
        character(*), intent(in) :: name
        integer, intent(inout) :: failures
        real(dp) :: observed, tolerance
        observed = sum(y)/real(size(y),dp)
        tolerance = 7.0_dp*sigma/sqrt(real(size(y),dp)) + &
            5.0e-4_dp*max(1.0_dp,abs(expected))
        if (abs(observed-expected) > tolerance) then
            print '(a,3es24.16)', trim(name)//': ', observed, expected, tolerance
            failures = failures + 1
        end if
    end subroutine check_mean

end program test_rng
