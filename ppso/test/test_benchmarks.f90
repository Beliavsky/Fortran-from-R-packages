program test_benchmarks
    use ppso, only : dp, rastrigin_function, ackley_function, griewank_function, sample_function
    implicit none
    real(dp) :: x(3)
    real(dp), parameter :: e = 2.7182818284590452353602874713526625_dp

    x = 0.0_dp
    call check(abs(rastrigin_function(x)+3.0_dp) < 1.0e-14_dp, "rastrigin zero")
    call check(abs(ackley_function(x)+20.0_dp+e) < 1.0e-13_dp, "ackley zero")
    call check(abs(griewank_function(x)) < 1.0e-14_dp, "griewank zero")
    call check(abs(sample_function(x)+6.0_dp) < 1.0e-14_dp, "sample zero")
    print *, "test_benchmarks: PASS"
contains
    subroutine check(ok, msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) error stop msg
    end subroutine check
end program test_benchmarks
