program ripple_example
    use subplex, only: dp, subplex_result, subplex_minimize
    implicit none

    type(subplex_result) :: result
    real(dp) :: x0(3)

    x0 = [0.1_dp, 3.0_dp, 2.0_dp]
    call subplex_minimize(ripple, x0, result, compute_hessian=.true.)
    write(*,'(a,3f14.7)') "x = ", result%x
    write(*,'(a,f14.8)') "f = ", result%value

contains

    function ripple(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f, r
        r = sqrt(sum(x*x))
        f = 1.0_dp - exp(-r*r)*cos(10.0_dp*r)**2
    end function ripple

end program ripple_example
