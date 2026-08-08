program test_all2one
    use soma, only : dp, soma_bounds, soma_result, bounds, soma_optimize, soma_set_seed
    implicit none

    type(soma_bounds) :: bnds
    type(soma_result) :: result
    integer :: i

    interface
        function sphere(x) result(value)
            use soma, only : dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function sphere
    end interface

    bnds = bounds([-5.0_dp, -5.0_dp], [5.0_dp, 5.0_dp])
    call soma_set_seed(1234)
    call soma_optimize(sphere, bnds, result)

    if (result%status /= 0) error stop 1
    if (result%cost(result%leader) > 0.1_dp) error stop 2
    if (maxval(abs(result%population(:,result%leader))) > 0.1_dp) error stop 3
    if (result%migrations /= 20) error stop 4
    if (size(result%history) /= result%migrations + 1) error stop 5
    if (size(result%evaluations) /= result%migrations + 1) error stop 6
    do i = 2, size(result%history)
        if (result%history(i) > result%history(i-1) + 1.0e-12_dp) error stop 7
    end do
end program test_all2one

function sphere(x) result(value)
    use soma, only : dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(x*x)
end function sphere
