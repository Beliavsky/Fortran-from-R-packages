program test_statistics
    use zigg, only : dp, zsetseed, zrnorm, zrexp, zrunif
    implicit none

    integer, parameter :: n = 100000
    real(dp), parameter :: tol = 5.0e-12_dp
    real(dp), allocatable :: x(:)
    real(dp) :: m, v

    call zsetseed(12345)

    x = zrnorm(n)
    m = sum(x) / real(n, dp)
    v = sum((x - m)**2) / real(n - 1, dp)
    call check(m, -0.0013092303662096288_dp, 'normal mean')
    call check(v, 1.0067727019522157_dp, 'normal variance')

    x = zrexp(n)
    m = sum(x) / real(n, dp)
    v = sum((x - m)**2) / real(n - 1, dp)
    call check(m, 1.0166745354931719_dp, 'exponential mean')
    call check(v, 1.1676628230956538_dp, 'exponential variance')
    if (minval(x) < 0.0_dp) error stop 'exponential support'

    x = zrunif(n)
    m = sum(x) / real(n, dp)
    v = sum((x - m)**2) / real(n - 1, dp)
    call check(m, 0.50115140948011661_dp, 'uniform mean')
    call check(v, 0.083646004972311297_dp, 'uniform variance')
    if (minval(x) < 0.0_dp .or. maxval(x) > 1.0_dp) error stop 'uniform support'

    print '(a)', 'test_statistics: PASS'

contains

    subroutine check(value, reference, label)
        real(dp), intent(in) :: value, reference
        character(len=*), intent(in) :: label
        if (abs(value - reference) > tol) then
            print '(a,2(1x,es24.16))', trim(label), value, reference
            error stop 'upstream stream-statistic mismatch'
        end if
    end subroutine check

end program test_statistics
