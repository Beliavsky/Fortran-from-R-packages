program test_compatibility
    use ppso, only : dp, dds_control, ppso_result, optim_dds
    implicit none
    type(dds_control) :: control
    type(ppso_result) :: legacy, corrected
    real(dp) :: lower(2), upper(2)

    lower = -1.0_dp
    upper = 1.0_dp
    control%number_of_particles = 1
    control%max_number_function_calls = 120
    control%seed = 77_8
    control%legacy_serial_prerun_omission = .true.
    call optim_dds(sphere, lower, upper, legacy, control)
    call check(legacy%actual_function_calls == legacy%function_calls+1, &
        "legacy DDS pre-run omission accounting")

    control%seed = 77_8
    control%legacy_serial_prerun_omission = .false.
    call optim_dds(sphere, lower, upper, corrected, control)
    call check(corrected%actual_function_calls == corrected%function_calls, &
        "corrected DDS call accounting")
    print *, "test_compatibility: PASS"
contains
    function sphere(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = sum(x*x)
    end function sphere
    subroutine check(ok, msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) error stop msg
    end subroutine check
end program test_compatibility
