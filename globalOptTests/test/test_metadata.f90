program test_metadata
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use global_opt_tests, only : dp, n_benchmarks, benchmark_names, get_default_bounds, get_problem_dimension, &
        get_global_opt, go_test
    implicit none
    real(dp), allocatable :: lower(:), upper(:)
    real(dp) :: value
    integer :: i, n, status

    if (n_benchmarks /= 50) error stop 'wrong benchmark count'
    do i = 1, n_benchmarks
        n = get_problem_dimension(trim(benchmark_names(i)), status)
        if (status /= 0 .or. n <= 0) error stop 'dimension lookup failed'
        call get_default_bounds(trim(benchmark_names(i)), lower, upper, status)
        if (status /= 0) error stop 'bounds lookup failed'
        if (size(lower) /= n .or. size(upper) /= n) error stop 'bounds size mismatch'
        if (any(upper < lower)) error stop 'invalid bounds'
        value = get_global_opt(trim(benchmark_names(i)), status)
        if (status /= 0 .or. ieee_is_nan(value)) error stop 'global optimum lookup failed'
        deallocate(lower, upper)
    end do

    value = go_test([0.0_dp], 'Ackleys', status=status)
    if (status /= 1 .or. .not. ieee_is_nan(value)) error stop 'dimension error path failed'
    value = go_test([0.0_dp], 'NotAFunction', status=status)
    if (status /= 2 .or. .not. ieee_is_nan(value)) error stop 'unknown-name path failed'

    print '(a)', 'test_metadata: PASS'
end program test_metadata
