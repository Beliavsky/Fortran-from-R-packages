program strategies_example
    use soma, only : dp, soma_bounds, soma_options, soma_result, bounds, all2one, t3a, pareto, &
                     soma_optimize, soma_set_seed
    implicit none

    type(soma_bounds) :: bnds
    type(soma_options) :: options
    type(soma_result) :: result

    interface
        function sphere(x) result(value)
            use soma, only : dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function sphere
    end interface

    bnds = bounds([-5.0_dp, -5.0_dp], [5.0_dp, 5.0_dp])

    call soma_set_seed(1234)
    options = all2one()
    call soma_optimize(sphere, bnds, result, options)
    print '(a,es12.4)', 'All To One: ', result%cost(result%leader)

    call soma_set_seed(1234)
    options = t3a()
    call soma_optimize(sphere, bnds, result, options)
    print '(a,es12.4)', 'T3A:       ', result%cost(result%leader)

    call soma_set_seed(1234)
    options = pareto()
    call soma_optimize(sphere, bnds, result, options)
    print '(a,es12.4)', 'Pareto:    ', result%cost(result%leader)
end program strategies_example

function sphere(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
