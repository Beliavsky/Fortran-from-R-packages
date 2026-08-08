program test_init_bookkeeping
    use soma, only : dp, soma_bounds, soma_options, soma_result, bounds, all2one, soma_optimize
    implicit none

    type(soma_bounds) :: bnds
    type(soma_options) :: options
    type(soma_result) :: result
    real(dp) :: init(2,3)

    interface
        function sphere(x) result(value)
            use soma, only : dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function sphere
    end interface

    init(:,1) = [2.0_dp, 2.0_dp]
    init(:,2) = [0.5_dp, -0.5_dp]
    init(:,3) = [1.0_dp, 0.0_dp]
    bnds = bounds([-5.0_dp, -5.0_dp], [5.0_dp, 5.0_dp])
    options = all2one(population_size=3, n_migrations=0)

    call soma_optimize(sphere, bnds, result, options, init)

    if (result%status /= 0) error stop 1
    if (result%migrations /= 0) error stop 2
    if (result%leader /= 2) error stop 3
    if (size(result%history) /= 1 .or. size(result%evaluations) /= 1) error stop 4
    if (result%evaluations(1) /= 0) error stop 5
    if (maxval(abs(result%population-init)) > 0.0_dp) error stop 6
end program test_init_bookkeeping

function sphere(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
