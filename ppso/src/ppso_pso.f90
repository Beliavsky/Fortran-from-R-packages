module ppso_pso
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use ppso_kinds, only : dp
    use ppso_types, only : pso_control, pso_state, ppso_result, objective_function, &
        ppso_monitor, PPSO_STOP_NONE, PPSO_STOP_MAX_CALLS, PPSO_STOP_MAX_ITER, &
        PPSO_STOP_CONVERGED, PPSO_STOP_INVALID_OBJECTIVE, stop_message
    use ppso_init, only : generate_initial_population, append_pending, validate_bounds
    implicit none
    private
    public :: optim_pso, init_pso_state, run_pso_state

contains

    subroutine init_pso_state(control, lower, upper, state, vmax, initial_estimates)
        type(pso_control), intent(in) :: control
        real(dp), intent(in) :: lower(:), upper(:)
        type(pso_state), intent(out) :: state
        real(dp), intent(in), optional :: vmax(:)
        real(dp), intent(in), optional :: initial_estimates(:,:)
        integer :: npar, npart

        call validate_bounds(lower, upper)
        npar = size(lower)
        npart = control%number_of_particles
        if (npart <= 0) error stop "optim_pso: number_of_particles must be positive"
        if (control%max_number_function_calls < npart) &
            error stop "optim_pso: max_number_function_calls must be at least number_of_particles"

        state%npar = npar
        state%npart = npart
        allocate(state%lower(npar), state%upper(npar), state%vmax(npar))
        state%lower = lower
        state%upper = upper
        if (present(vmax)) then
            if (size(vmax) /= npar) error stop "optim_pso: Vmax dimension mismatch"
            if (any(vmax < 0.0_dp)) error stop "optim_pso: Vmax must be nonnegative"
            state%vmax = vmax
        else
            state%vmax = (upper-lower)/3.0_dp
        end if

        call state%rng%seed(control%seed)
        if (present(initial_estimates)) then
            call generate_initial_population(state%rng, lower, upper, npart, control%lhs_init, &
                initial_estimates, state%x, state%pending_initial)
        else
            call generate_initial_population(state%rng, lower, upper, npart, control%lhs_init, &
                x=state%x, pending=state%pending_initial)
        end if

        allocate(state%v(npar,npart), state%fitness_x(npart))
        allocate(state%x_lbest(npar,npart), state%fitness_lbest(npart))
        allocate(state%x_gbest(npar), state%function_calls(npart))
        allocate(state%fixed_until_evaluated(npart))
        state%v = 0.0_dp
        state%fitness_x = huge(1.0_dp)
        state%x_lbest = state%x
        state%fitness_lbest = huge(1.0_dp)
        state%x_gbest = huge(1.0_dp)
        state%function_calls = 0
        state%fixed_until_evaluated = .true.
        state%fitness_gbest = huge(1.0_dp)
        state%fitness_itbest = huge(1.0_dp)
        state%it_last_improvement = 0
        state%stop_code = PPSO_STOP_NONE
    end subroutine init_pso_state

    subroutine optim_pso(objective, lower, upper, result, control, vmax, initial_estimates, monitor, final_state)
        procedure(objective_function) :: objective
        real(dp), intent(in) :: lower(:), upper(:)
        type(ppso_result), intent(out) :: result
        type(pso_control), intent(in), optional :: control
        real(dp), intent(in), optional :: vmax(:)
        real(dp), intent(in), optional :: initial_estimates(:,:)
        procedure(ppso_monitor), optional :: monitor
        type(pso_state), intent(out), optional :: final_state
        type(pso_control) :: ctl
        type(pso_state) :: state

        if (present(control)) ctl = control
        if (present(vmax)) then
            if (present(initial_estimates)) then
                call init_pso_state(ctl, lower, upper, state, vmax, initial_estimates)
            else
                call init_pso_state(ctl, lower, upper, state, vmax=vmax)
            end if
        else
            if (present(initial_estimates)) then
                call init_pso_state(ctl, lower, upper, state, initial_estimates=initial_estimates)
            else
                call init_pso_state(ctl, lower, upper, state)
            end if
        end if
        call run_pso_state(objective, ctl, state, result, monitor)
        if (present(final_state)) final_state = state
    end subroutine optim_pso

    subroutine run_pso_state(objective, control, state, result, monitor)
        procedure(objective_function) :: objective
        type(pso_control), intent(in) :: control
        type(pso_state), intent(inout) :: state
        type(ppso_result), intent(out) :: result
        procedure(ppso_monitor), optional :: monitor
        integer :: i
        logical, allocatable :: completed(:)

        allocate(completed(state%npart))

        do while (state%stop_code == PPSO_STOP_NONE)
            if (control%wait_complete_iteration) then
                completed = .false.
                do i = 1, state%npart
                    ! The R wait_complete_iteration path evaluates the full swarm once
                    ! a generation has started, even if this overshoots the call budget.
                    state%fitness_x(i) = objective(state%x(:,i))
                    state%function_calls(i) = state%function_calls(i) + 1
                    state%fixed_until_evaluated(i) = .false.
                    completed(i) = .true.
                    if (ieee_is_nan(state%fitness_x(i))) then
                        state%stop_code = PPSO_STOP_INVALID_OBJECTIVE
                        exit
                    end if
                end do
                if (any(completed)) call update_pso(control, state, completed)
            else
                do i = 1, state%npart
                    if (state%stop_code /= PPSO_STOP_NONE) exit
                    if (sum(state%function_calls) >= control%max_number_function_calls) then
                        state%stop_code = PPSO_STOP_MAX_CALLS
                        exit
                    end if
                    state%fitness_x(i) = objective(state%x(:,i))
                    state%function_calls(i) = state%function_calls(i) + 1
                    state%fixed_until_evaluated(i) = .false.
                    if (ieee_is_nan(state%fitness_x(i))) then
                        state%stop_code = PPSO_STOP_INVALID_OBJECTIVE
                        exit
                    end if
                    completed = .false.
                    completed(i) = .true.
                    call update_pso(control, state, completed)
                end do
            end if

            if (present(monitor)) then
                call monitor(minval(state%function_calls), sum(state%function_calls), &
                    state%fitness_gbest, state%x_gbest)
            end if
        end do

        result%value = state%fitness_gbest
        allocate(result%par(state%npar))
        result%par = state%x_gbest
        result%function_calls = sum(state%function_calls)
        result%actual_function_calls = result%function_calls
        result%iterations = minval(state%function_calls)
        result%stop_code = state%stop_code
        result%break_flag = stop_message(state%stop_code)
    end subroutine run_pso_state

    subroutine update_pso(control, state, completed)
        type(pso_control), intent(in) :: control
        type(pso_state), intent(inout) :: state
        logical, intent(inout) :: completed(:)
        logical, allocatable :: improved(:)
        integer :: i, j, min_index
        real(dp) :: min_fitness, scale_factor, denom, rel_improvement
        real(dp), allocatable :: r1(:), r2(:)
        logical :: used_pending

        if (size(completed) /= state%npart) error stop "update_pso: completed size mismatch"
        allocate(improved(state%npart), r1(state%npar), r2(state%npar))
        improved = state%fitness_x < state%fitness_lbest

        do i = 1, state%npart
            if (completed(i) .and. improved(i)) then
                state%fitness_lbest(i) = state%fitness_x(i)
                state%x_lbest(:,i) = state%x(:,i)
            end if
        end do

        min_fitness = huge(1.0_dp)
        min_index = 0
        do i = 1, state%npart
            if (completed(i) .and. state%fitness_x(i) < min_fitness) then
                min_fitness = state%fitness_x(i)
                min_index = i
            end if
        end do
        if (min_index > 0 .and. min_fitness < state%fitness_gbest) then
            state%fitness_gbest = min_fitness
            state%x_gbest = state%x(:,min_index)
            do i = 1, state%npart
                if (.not. state%fixed_until_evaluated(i)) completed(i) = .true.
            end do
        end if

        denom = max(1.0e-5_dp, abs(state%fitness_gbest))
        rel_improvement = abs((state%fitness_itbest-state%fitness_gbest)/denom)
        if ((state%fitness_itbest-state%fitness_gbest > control%abstol) .and. &
            (rel_improvement > control%reltol)) then
            state%it_last_improvement = maxval(state%function_calls)
            state%fitness_itbest = state%fitness_gbest
        else if (minval(state%function_calls)-state%it_last_improvement >= control%max_wait_iterations) then
            state%stop_code = PPSO_STOP_CONVERGED
        end if

        if (minval(state%function_calls) >= control%max_number_of_iterations) &
            state%stop_code = PPSO_STOP_MAX_ITER
        if (sum(state%function_calls) >= control%max_number_function_calls) &
            state%stop_code = PPSO_STOP_MAX_CALLS
        if (state%stop_code /= PPSO_STOP_NONE) return

        do i = 1, state%npart
            if (.not. completed(i)) cycle
            call append_pending(state%pending_initial, state%x(:,i), used_pending)
            if (used_pending) then
                state%fixed_until_evaluated(i) = .true.
                state%v(:,i) = 0.0_dp
            else
                do j = 1, state%npar
                    r1(j) = state%rng%uniform()
                    r2(j) = state%rng%uniform()
                end do
                state%v(:,i) = control%w*state%v(:,i) &
                    + control%c1*r1*(state%x_lbest(:,i)-state%x(:,i)) &
                    + control%c2*r2*(state%x_gbest-state%x(:,i))
                scale_factor = 1.0_dp
                do j = 1, state%npar
                    if (abs(state%v(j,i)) > 0.0_dp) then
                        scale_factor = min(scale_factor, abs(state%vmax(j)/state%v(j,i)))
                    end if
                end do
                state%v(:,i) = state%v(:,i)*min(1.0_dp,scale_factor)
                state%x(:,i) = state%x(:,i) + state%v(:,i)
                state%fixed_until_evaluated(i) = .false.
            end if
            state%x(:,i) = max(state%x(:,i), state%lower)
            state%x(:,i) = min(state%x(:,i), state%upper)
            state%fitness_x(i) = huge(1.0_dp)
        end do
    end subroutine update_pso

end module ppso_pso
