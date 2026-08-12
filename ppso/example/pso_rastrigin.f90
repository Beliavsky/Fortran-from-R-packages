program pso_rastrigin
    use ppso, only : dp, pso_control, ppso_result, optim_pso, rastrigin_function
    implicit none
    type(pso_control) :: control
    type(ppso_result) :: result
    real(dp) :: lower(2), upper(2)

    lower = -1.0_dp
    upper = 1.0_dp
    control%number_of_particles = 40
    control%max_number_of_iterations = 80
    control%max_number_function_calls = 3200
    control%seed = 20260810_8

    call optim_pso(rastrigin_function, lower, upper, result, control)
    write(*,'(a,es16.8)') "value = ", result%value
    write(*,'(a,*(f12.7,1x))') "par   = ", result%par
    write(*,'(a,i0)') "calls = ", result%function_calls
    write(*,'(a,a)') "stop  = ", trim(result%break_flag)
end program pso_rastrigin
