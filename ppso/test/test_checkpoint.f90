program test_checkpoint
    use ppso, only : dp, pso_control, pso_state, ppso_result, init_pso_state, run_pso_state, &
        save_pso_state, load_pso_state
    implicit none
    type(pso_control) :: control
    type(pso_state) :: a, b
    type(ppso_result) :: ra, rb
    real(dp) :: lower(3), upper(3)
    character(len=*), parameter :: filename = "ppso_checkpoint_test.bin"

    lower = -2.0_dp
    upper = 2.0_dp
    control%number_of_particles = 20
    control%max_number_of_iterations = 30
    control%max_number_function_calls = 600
    control%seed = 12345_8

    call init_pso_state(control, lower, upper, a)
    call save_pso_state(filename, a)
    call load_pso_state(filename, b)
    call run_pso_state(sphere, control, a, ra)
    call run_pso_state(sphere, control, b, rb)
    call check(abs(ra%value-rb%value) <= 0.0_dp, "checkpoint value mismatch")
    call check(all(abs(ra%par-rb%par) <= 0.0_dp), "checkpoint parameter mismatch")
    call delete_file(filename)
    print *, "test_checkpoint: PASS"
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

    subroutine delete_file(name)
        character(len=*), intent(in) :: name
        integer :: u, ios
        open(newunit=u, file=name, status="old", action="readwrite", iostat=ios)
        if (ios == 0) close(u, status="delete")
    end subroutine delete_file
end program test_checkpoint
