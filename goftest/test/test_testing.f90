module test_cdf_module
    use goftest, only : dp
    implicit none
contains
    function normal_cdf(x) result(p)
        real(dp), intent(in) :: x
        real(dp) :: p
        p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
    end function normal_cdf
end module test_cdf_module

program test_testing
    use goftest, only : dp, gof_result, ad_test, cvm_test, ad_test_values, cvm_test_values
    use test_cdf_module, only : normal_cdf
    implicit none

    real(dp) :: x(8), u(12)
    integer :: g(12), i
    type(gof_result) :: a, c, ab, cb

    x = [-1.2815515655446004_dp, -0.8416212335729143_dp, -0.5244005127080409_dp, &
        -0.2533471031357997_dp, 0.0_dp, 0.2533471031357997_dp, &
        0.8416212335729143_dp, 1.2815515655446004_dp]
    a = ad_test(x, normal_cdf)
    c = cvm_test(x, normal_cdf)
    if (a%status /= 0 .or. c%status /= 0) error stop 1
    if (a%p_value <= 0.0_dp .or. a%p_value > 1.0_dp) error stop 1
    if (c%p_value <= 0.0_dp .or. c%p_value > 1.0_dp) error stop 1

    do i = 1, 12
        u(i) = (real(i, dp) - 0.35_dp) / 12.3_dp
        g(i) = 1 + mod(i - 1, 3)
    end do
    ab = ad_test_values(u, estimated=.true., groups=g)
    cb = cvm_test_values(u, estimated=.true., groups=g)
    if (ab%status /= 0 .or. cb%status /= 0) error stop 1
    if (ab%groups /= 3 .or. cb%groups /= 3) error stop 1
    if (ab%p_value < 0.0_dp .or. ab%p_value > 1.0_dp) error stop 1
    if (cb%p_value < 0.0_dp .or. cb%p_value > 1.0_dp) error stop 1

    print '(a)', 'test_testing: PASS'
end program test_testing
