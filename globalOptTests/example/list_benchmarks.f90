program list_benchmarks
    use global_opt_tests, only : n_benchmarks, benchmark_names, get_problem_dimension, get_global_opt
    implicit none
    integer :: i

    print '(a)', 'Name             dim   published optimum'
    do i = 1, n_benchmarks
        print '(a16,1x,i3,1x,es16.8)', benchmark_names(i), &
            get_problem_dimension(trim(benchmark_names(i))), get_global_opt(trim(benchmark_names(i)))
    end do
end program list_benchmarks
