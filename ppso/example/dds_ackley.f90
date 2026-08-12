program dds_ackley
    use ppso, only : dp, dds_control, ppso_result, optim_dds, ackley_function
    implicit none
    type(dds_control) :: control
    type(ppso_result) :: result
    real(dp) :: lower(4), upper(4)

    lower = -5.0_dp
    upper = 5.0_dp
    control%number_of_particles = 2
    control%max_number_function_calls = 1600
    control%part_xchange = 2
    control%legacy_serial_prerun_omission = .false.
    control%seed = 31415926_8

    call optim_dds(ackley_function, lower, upper, result, control)
    write(*,'(a,es16.8)') "value = ", result%value
    write(*,'(a,*(f12.7,1x))') "par   = ", result%par
    write(*,'(a,i0)') "calls = ", result%function_calls
    write(*,'(a,a)') "stop  = ", trim(result%break_flag)
end program dds_ackley
