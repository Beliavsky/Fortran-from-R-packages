program integer_map_example
    use deoptim, only : dp, i8, de_control, de_result, deoptim_solve
    implicit none
    type(de_control) :: control
    type(de_result) :: result
    real(dp) :: lower(3), upper(3)

    lower = -10.0_dp
    upper = 10.0_dp
    control = de_control()
    control%np = 50
    control%itermax = 150
    control%seed = 123_i8

    call deoptim_solve(objective, lower, upper, result, control, map=integer_map)
    write(*,'(a,es16.8)') "best value: ", result%bestval
    write(*,'(a,*(1x,f8.2))') "integer member:", result%bestmem
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (x(1)-2.0_dp)**2 + (x(2)+3.0_dp)**2 + (x(3)-5.0_dp)**2
    end function objective

    subroutine integer_map(x)
        real(dp), intent(inout) :: x(:)
        x = real(nint(x), dp)
    end subroutine integer_map
end program integer_map_example
