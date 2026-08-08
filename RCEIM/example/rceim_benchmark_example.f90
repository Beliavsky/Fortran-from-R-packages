program rceim_benchmark_example
    use rceim, only : dp, rceim_options, rceim_result, ceim_optimize
    use rceim_benchmarks, only : test_fun_optimization_2d
    implicit none
    type(rceim_options) :: opt
    type(rceim_result) :: res
    real(dp) :: lower(2), upper(2)

    lower = [-10.0_dp,-10.0_dp]
    upper = [ 10.0_dp, 10.0_dp]
    opt%n_total = 700
    opt%n_elite = 175
    opt%n_super = 2
    opt%epsilon = 0.03_dp
    opt%max_iter = 40
    opt%seed = 42
    call ceim_optimize(test_fun_optimization_2d, lower, upper, res, opt)
    write(*,'(a,2f12.6)') 'best x = ', res%x
    write(*,'(a,es14.6)') 'best value = ', res%value
end program rceim_benchmark_example
