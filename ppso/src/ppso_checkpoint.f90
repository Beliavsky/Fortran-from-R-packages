module ppso_checkpoint
    use ppso_kinds, only : dp
    use ppso_types, only : pso_state, dds_state
    implicit none
    private
    public :: save_pso_state, load_pso_state, save_dds_state, load_dds_state

contains

    subroutine save_pso_state(filename, state)
        character(len=*), intent(in) :: filename
        type(pso_state), intent(in) :: state
        integer :: u, npending
        character(len=16), parameter :: magic = "PPSO_PSO_V1     "

        npending = size(state%pending_initial,2)
        open(newunit=u, file=filename, access="stream", form="unformatted", status="replace", action="write")
        write(u) magic
        write(u) state%npar, state%npart, npending
        write(u) state%fitness_gbest, state%fitness_itbest, state%it_last_improvement, state%stop_code
        write(u) state%rng%state, state%rng%have_spare, state%rng%spare
        write(u) state%lower, state%upper, state%vmax
        write(u) state%x, state%v, state%fitness_x
        write(u) state%x_lbest, state%fitness_lbest, state%x_gbest
        if (npending > 0) write(u) state%pending_initial
        write(u) state%function_calls, state%fixed_until_evaluated
        close(u)
    end subroutine save_pso_state

    subroutine load_pso_state(filename, state)
        character(len=*), intent(in) :: filename
        type(pso_state), intent(out) :: state
        integer :: u, npending
        character(len=16) :: magic

        open(newunit=u, file=filename, access="stream", form="unformatted", status="old", action="read")
        read(u) magic
        if (magic /= "PPSO_PSO_V1     ") error stop "load_pso_state: incompatible checkpoint"
        read(u) state%npar, state%npart, npending
        allocate(state%lower(state%npar), state%upper(state%npar), state%vmax(state%npar))
        allocate(state%x(state%npar,state%npart), state%v(state%npar,state%npart), state%fitness_x(state%npart))
        allocate(state%x_lbest(state%npar,state%npart), state%fitness_lbest(state%npart))
        allocate(state%x_gbest(state%npar), state%pending_initial(state%npar,npending))
        allocate(state%function_calls(state%npart), state%fixed_until_evaluated(state%npart))
        read(u) state%fitness_gbest, state%fitness_itbest, state%it_last_improvement, state%stop_code
        read(u) state%rng%state, state%rng%have_spare, state%rng%spare
        read(u) state%lower, state%upper, state%vmax
        read(u) state%x, state%v, state%fitness_x
        read(u) state%x_lbest, state%fitness_lbest, state%x_gbest
        if (npending > 0) read(u) state%pending_initial
        read(u) state%function_calls, state%fixed_until_evaluated
        close(u)
    end subroutine load_pso_state

    subroutine save_dds_state(filename, state)
        character(len=*), intent(in) :: filename
        type(dds_state), intent(in) :: state
        integer :: u, npending
        character(len=16), parameter :: magic = "PPSO_DDS_V1     "

        npending = size(state%pending_initial,2)
        open(newunit=u, file=filename, access="stream", form="unformatted", status="replace", action="write")
        write(u) magic
        write(u) state%npar, state%npart, npending
        write(u) state%fitness_gbest, state%fitness_itbest, state%it_last_improvement
        write(u) state%init_function_calls, state%actual_init_calls, state%main_max_calls, state%stop_code
        write(u) state%rng%state, state%rng%have_spare, state%rng%spare
        write(u) state%lower, state%upper
        write(u) state%x, state%v, state%fitness_x
        write(u) state%x_lbest, state%fitness_lbest, state%x_gbest
        if (npending > 0) write(u) state%pending_initial
        write(u) state%function_calls, state%futile_iter_count
        close(u)
    end subroutine save_dds_state

    subroutine load_dds_state(filename, state)
        character(len=*), intent(in) :: filename
        type(dds_state), intent(out) :: state
        integer :: u, npending
        character(len=16) :: magic

        open(newunit=u, file=filename, access="stream", form="unformatted", status="old", action="read")
        read(u) magic
        if (magic /= "PPSO_DDS_V1     ") error stop "load_dds_state: incompatible checkpoint"
        read(u) state%npar, state%npart, npending
        allocate(state%lower(state%npar), state%upper(state%npar))
        allocate(state%x(state%npar,state%npart), state%v(state%npar,state%npart), state%fitness_x(state%npart))
        allocate(state%x_lbest(state%npar,state%npart), state%fitness_lbest(state%npart))
        allocate(state%x_gbest(state%npar), state%pending_initial(state%npar,npending))
        allocate(state%function_calls(state%npart), state%futile_iter_count(state%npart))
        read(u) state%fitness_gbest, state%fitness_itbest, state%it_last_improvement
        read(u) state%init_function_calls, state%actual_init_calls, state%main_max_calls, state%stop_code
        read(u) state%rng%state, state%rng%have_spare, state%rng%spare
        read(u) state%lower, state%upper
        read(u) state%x, state%v, state%fitness_x
        read(u) state%x_lbest, state%fitness_lbest, state%x_gbest
        if (npending > 0) read(u) state%pending_initial
        read(u) state%function_calls, state%futile_iter_count
        close(u)
    end subroutine load_dds_state

end module ppso_checkpoint
