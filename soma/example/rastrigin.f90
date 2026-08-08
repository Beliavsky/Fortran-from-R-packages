program rastrigin_example
    use soma, only : dp, soma_bounds, soma_result, bounds, soma_optimize, soma_set_seed
    implicit none

    type(soma_bounds) :: bnds
    type(soma_result) :: result

    interface
        function rastrigin(x) result(value)
            use soma, only : dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function rastrigin
    end interface

    bnds = bounds([-5.12_dp, -5.12_dp], [5.12_dp, 5.12_dp])
    call soma_set_seed(2026)
    call soma_optimize(rastrigin, bnds, result)

    print '(a,es14.6)', 'best cost: ', result%cost(result%leader)
    print '(a,*(f12.6,1x))', 'best point: ', result%population(:,result%leader)
    print '(a,i0)', 'migrations: ', result%migrations
end program rastrigin_example

function rastrigin(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp) :: pi_value

    pi_value = acos(-1.0_dp)
    value = 10.0_dp * real(size(x),dp) + sum(x*x - 10.0_dp*cos(2.0_dp*pi_value*x))
end function rastrigin
