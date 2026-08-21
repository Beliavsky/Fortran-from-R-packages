program rosenbrock_example
    use deoptim, only : dp, i8, de_control, de_result, deoptim_solve
    implicit none
    type(de_control) :: control
    type(de_result) :: result
    real(dp) :: lower(2), upper(2)

    lower = -10.0_dp
    upper = 10.0_dp
    control = de_control()
    control%np = 50
    control%itermax = 400
    control%strategy = 2
    control%cr = 0.9_dp
    control%seed = 21_i8
    control%trace = 50

    call deoptim_solve(rosenbrock, lower, upper, result, control)
    write(*,'(a,es16.8)') "best value: ", result%bestval
    write(*,'(a,*(1x,es16.8))') "best member:", result%bestmem
    write(*,'(a,i0)') "generations: ", result%iter
    write(*,'(a,i0)') "function evaluations: ", result%nfeval
contains
    function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 100.0_dp*(x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2
    end function rosenbrock
end program rosenbrock_example
