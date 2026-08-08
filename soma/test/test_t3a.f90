program test_t3a
    use soma, only : dp, soma_bounds, soma_options, soma_result, bounds, t3a, soma_optimize, soma_set_seed
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
    options = t3a()
    call soma_set_seed(1234)
    call soma_optimize(sphere, bnds, result, options)

    if (result%status /= 0) error stop 1
    if (result%cost(result%leader) > 1.0e-3_dp) error stop 2
    if (result%migrations /= options%n_migrations) error stop 3
    if (any(result%population < -5.0_dp) .or. any(result%population > 5.0_dp)) error stop 4
end program test_t3a

function sphere(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
