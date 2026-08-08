program test_validation
    use soma, only : dp, soma_bounds, soma_options, soma_result, bounds, t3a, soma_optimize
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

    bnds = bounds([-1.0_dp], [1.0_dp])
    options = t3a(population_size=5, leader_pool_size=6)
    call soma_optimize(sphere, bnds, result, options)
    if (result%status == 0) error stop 1
end program test_validation

function sphere(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
