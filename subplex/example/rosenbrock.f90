program rosenbrock_example
    use subplex, only: dp, subplex_result, subplex_minimize
    implicit none

    type(subplex_result) :: result
    real(dp) :: x0(2)

    x0 = [11.0_dp, -33.0_dp]
    call subplex_minimize(rosen, x0, result)
    write(*,'(a,2f16.8)') "x = ", result%x
    write(*,'(a,es16.8)') "f = ", result%value
    write(*,'(a,i0)') "evaluations = ", result%counts
    write(*,'(a,i0,2a)') "convergence = ", result%convergence, ": ", &
        result%message

contains

    function rosen(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
    end function rosen

end program rosenbrock_example
