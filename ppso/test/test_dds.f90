program test_dds
    use ppso, only : dp, dds_control, ppso_result, optim_dds
    implicit none
    type(dds_control) :: control
    type(ppso_result) :: result
    real(dp) :: lower(4), upper(4)
    integer :: mode

    lower = -4.0_dp
    upper = 4.0_dp
    control%number_of_particles = 4
    control%max_number_function_calls = 1400
    control%legacy_serial_prerun_omission = .false.

    do mode = 0, 3
        control%part_xchange = mode
        control%seed = 1000_8 + int(mode,kind=8)
        call optim_dds(sphere, lower, upper, result, control)
        call check(result%value < 2.0e-3_dp, "DDS mode failed on sphere")
        call check(all(result%par >= lower .and. result%par <= upper), "DDS bounds")
    end do
    print *, "test_dds: PASS"
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
end program test_dds
