program basic_usage
    use global_opt_tests, only : dp, go_test, get_default_bounds, get_global_opt
    implicit none
    real(dp), allocatable :: lower(:), upper(:), x(:)
    real(dp) :: value

    call get_default_bounds('Branin', lower, upper)
    allocate(x(size(lower)))
    x = [-3.14159265359_dp, 12.275_dp]
    value = go_test(x, 'Branin')

    print '(a,2(1x,f12.6))', 'Branin point:', x
    print '(a,1x,es16.8)', 'value:', value
    print '(a,1x,es16.8)', 'published global optimum:', get_global_opt('Branin')
end program basic_usage
