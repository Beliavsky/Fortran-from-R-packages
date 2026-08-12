program test_pso
    use ppso, only : dp, pso_control, ppso_result, optim_pso, rastrigin_function
    implicit none
    type(pso_control) :: control
    type(ppso_result) :: a, b
    real(dp) :: lower(2), upper(2)

    lower = -1.0_dp
    upper = 1.0_dp
    control%number_of_particles = 40
    control%max_number_of_iterations = 100
    control%max_number_function_calls = 4000
    control%seed = 9918273_8

    control%wait_complete_iteration = .false.
    call optim_pso(rastrigin_function, lower, upper, a, control)
    call check(a%value < -1.9999_dp, "asynchronous PSO did not reach Rastrigin optimum")
    call check(all(a%par >= lower .and. a%par <= upper), "asynchronous PSO bounds")

    control%wait_complete_iteration = .true.
    control%seed = 9918273_8
    call optim_pso(rastrigin_function, lower, upper, b, control)
    call check(b%value < -1.999_dp, "synchronous PSO did not reach Rastrigin optimum")
    call check(b%actual_function_calls == b%function_calls, "PSO call accounting")
    print *, "test_pso: PASS"
contains
    subroutine check(ok, msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) error stop msg
    end subroutine check
end program test_pso
