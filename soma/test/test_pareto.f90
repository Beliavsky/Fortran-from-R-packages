program test_pareto
    use soma, only : dp, soma_bounds, soma_options, soma_result, bounds, pareto, soma_optimize, soma_set_seed
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
    options = pareto()
    call soma_set_seed(1234)
    call soma_optimize(sphere, bnds, result, options)

    if (result%status /= 0) error stop 1
    if (result%cost(result%leader) > 0.1_dp) error stop 2
    if (result%history(size(result%history)) > result%history(1)) error stop 3
    if (any(result%population < -5.0_dp) .or. any(result%population > 5.0_dp)) error stop 4
end program test_pareto

function sphere(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
