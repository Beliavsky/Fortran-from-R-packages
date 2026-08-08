program test_benchmarks
    use rceim_kinds, only : dp
    use rceim_benchmarks, only : test_fun_optimization, test_fun_optimization_2d
    implicit none
    real(dp) :: x1(1), x2(2), expected
    x1 = [0.0_dp]
    expected = exp(-4.0_dp) + 0.9_dp*exp(-4.0_dp) + 0.25_dp
    if (abs(test_fun_optimization(x1)-expected) > 1.0e-13_dp) error stop 1
    x2 = [0.0_dp,0.0_dp]
    expected = 20.0_dp/50.0_dp - 20.0_dp/90.0_dp - exp(-4.0_dp) &
             - 0.9_dp*exp(-4.0_dp) - 0.25_dp + 0.5_dp
    if (abs(test_fun_optimization_2d(x2)-expected) > 1.0e-13_dp) error stop 2
end program test_benchmarks
