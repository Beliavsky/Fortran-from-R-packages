module ppso_dds
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use ppso_kinds, only : dp
    use ppso_types, only : dds_control, dds_state, ppso_result, objective_function, &
        ppso_monitor, PPSO_STOP_NONE, PPSO_STOP_MAX_CALLS, PPSO_STOP_CONVERGED, &
        PPSO_STOP_INVALID_OBJECTIVE, stop_message
    use ppso_init, only : generate_initial_population, append_pending, validate_bounds
    implicit none
    private
    public :: optim_dds, init_dds_state, run_dds_state

contains

    subroutine optim_dds(objective, lower, upper, result, control, initial_estimates, monitor, final_state)
        procedure(objective_function) :: objective
        real(dp), intent(in) :: lower(:), upper(:)
        type(ppso_result), intent(out) :: result
        type(dds_control), intent(in), optional :: control
        real(dp), intent(in), optional :: initial_estimates(:,:)
        procedure(ppso_monitor), optional :: monitor
        type(dds_state), intent(out), optional :: final_state
        type(dds_control) :: ctl
        type(dds_state) :: state

        if (present(control)) ctl = control
        if (present(initial_estimates)) then
            call init_dds_state(objective, ctl, lower, upper, state, initial_estimates)
        else
            call init_dds_state(objective, ctl, lower, upper, state)
        end if
        if (state%stop_code == PPSO_STOP_NONE) call run_dds_state(objective, ctl, state, result, monitor)
        if (state%stop_code /= PPSO_STOP_NONE .and. .not. allocated(result%par)) then
            result%value = state%fitness_gbest
            allocate(result%par(state%npar))
            result%par = state%x_gbest
            result%function_calls = state%init_function_calls + sum(state%function_calls)
            result%actual_function_calls = state%actual_init_calls + sum(state%function_calls)
            result%iterations = minval(state%function_calls)
            result%stop_code = state%stop_code
            result%break_flag = stop_message(state%stop_code)
        end if
        if (present(final_state)) final_state = state
    end subroutine optim_dds

    subroutine init_dds_state(objective, control, lower, upper, state, initial_estimates)
        procedure(objective_function) :: objective
        type(dds_control), intent(in) :: control
        real(dp), intent(in) :: lower(:), upper(:)
        type(dds_state), intent(out) :: state
        real(dp), intent(in), optional :: initial_estimates(:,:)
        integer :: npar, npart, init_calls, npre, i, j, nkeep, offset
        integer, allocatable :: order(:)
        real(dp), allocatable :: prex(:,:), pending(:,:), fpre(:)
        logical :: legacy_drop

        call validate_bounds(lower, upper)
        npar = size(lower)
        npart = control%number_of_particles
        if (npart <= 0) error stop "optim_dds: number_of_particles must be positive"
        if (control%max_number_function_calls < npart) &
            error stop "optim_dds: max_number_function_calls must be at least number_of_particles"
        if (control%part_xchange < 0 .or. control%part_xchange > 3) &
            error stop "optim_dds: part_xchange must be 0, 1, 2, or 3"

        state%npar = npar
        state%npart = npart
        allocate(state%lower(npar), state%upper(npar))
        state%lower = lower
        state%upper = upper
        call state%rng%seed(control%seed)

        init_calls = ceiling(max(0.005_dp*real(control%max_number_function_calls,dp), 5.0_dp))
        npre = max(npart, init_calls)
        if (present(initial_estimates)) then
            call generate_initial_population(state%rng, lower, upper, npre, control%lhs_init, &
                initial_estimates, prex, pending)
        else
            call generate_initial_population(state%rng, lower, upper, npre, control%lhs_init, &
                x=prex, pending=pending)
        end if

        allocate(fpre(npre), order(npre))
        do i = 1, npre
            fpre(i) = objective(prex(:,i))
            state%actual_init_calls = state%actual_init_calls + 1
            if (ieee_is_nan(fpre(i))) then
                state%stop_code = PPSO_STOP_INVALID_OBJECTIVE
                exit
            end if
        end do
        if (state%stop_code /= PPSO_STOP_NONE) then
            allocate(state%x_gbest(npar))
            state%x_gbest = huge(1.0_dp)
            return
        end if

        legacy_drop = control%legacy_serial_prerun_omission .and. (npre-1 >= npart)
        offset = merge(1, 0, legacy_drop)
        state%init_function_calls = npre-offset
        state%main_max_calls = control%max_number_function_calls-state%init_function_calls
        if (state%main_max_calls <= 0) then
            state%stop_code = PPSO_STOP_MAX_CALLS
            allocate(state%x_gbest(npar))
            state%x_gbest = huge(1.0_dp)
            return
        end if

        do i = 1, npre
            order(i) = i
        end do
        call sort_indices_by_value(fpre, order, 1+offset, npre)
        nkeep = npart

        allocate(state%x(npar,npart), state%v(npar,npart), state%fitness_x(npart))
        allocate(state%x_lbest(npar,npart), state%fitness_lbest(npart))
        allocate(state%x_gbest(npar), state%function_calls(npart), state%futile_iter_count(npart))
        do i = 1, nkeep
            j = order(offset+i)
            state%x_lbest(:,i) = prex(:,j)
            state%fitness_lbest(i) = fpre(j)
        end do
        state%x = state%x_lbest
        state%v = 0.0_dp
        state%fitness_x = state%fitness_lbest
        state%function_calls = 0
        state%futile_iter_count = 0
        state%pending_initial = pending
        state%fitness_gbest = minval(state%fitness_lbest)
        j = minloc(state%fitness_lbest, dim=1)
        state%x_gbest = state%x_lbest(:,j)
        state%fitness_itbest = state%fitness_gbest
        state%it_last_improvement = 0

        ! optim_dds calls update_tasklist_dds(loop_counter=0) before the main search.
        call propose_dds(control, state)
    end subroutine init_dds_state

    subroutine run_dds_state(objective, control, state, result, monitor)
        procedure(objective_function) :: objective
        type(dds_control), intent(in) :: control
        type(dds_state), intent(inout) :: state
        type(ppso_result), intent(out) :: result
        procedure(ppso_monitor), optional :: monitor
        integer :: i

        do while (state%stop_code == PPSO_STOP_NONE)
            do i = 1, state%npart
                state%fitness_x(i) = objective(state%x(:,i))
                state%function_calls(i) = state%function_calls(i)+1
                if (ieee_is_nan(state%fitness_x(i))) then
                    state%stop_code = PPSO_STOP_INVALID_OBJECTIVE
                    exit
                end if
            end do
            if (state%stop_code /= PPSO_STOP_NONE) exit
            call update_dds(control, state)
            if (present(monitor)) then
                call monitor(minval(state%function_calls), &
                    state%init_function_calls+sum(state%function_calls), &
                    state%fitness_gbest, state%x_gbest)
            end if
        end do

        result%value = state%fitness_gbest
        allocate(result%par(state%npar))
        result%par = state%x_gbest
        result%function_calls = state%init_function_calls + sum(state%function_calls)
        result%actual_function_calls = state%actual_init_calls + sum(state%function_calls)
        result%iterations = minval(state%function_calls)
        result%stop_code = state%stop_code
        result%break_flag = stop_message(state%stop_code)
    end subroutine run_dds_state

    subroutine update_dds(control, state)
        type(dds_control), intent(in) :: control
        type(dds_state), intent(inout) :: state
        logical, allocatable :: improved(:), relocate(:)
        integer :: i, j, idx, recent
        real(dp) :: min_fitness, denom, rel_improvement

        allocate(improved(state%npart), relocate(state%npart))
        improved = state%fitness_x < state%fitness_lbest
        do i = 1, state%npart
            if (improved(i)) then
                state%futile_iter_count(i) = 0
                state%fitness_lbest(i) = state%fitness_x(i)
                state%x_lbest(:,i) = state%x(:,i)
            else
                state%futile_iter_count(i) = state%futile_iter_count(i)+1
            end if
        end do

        if (any(improved)) then
            min_fitness = minval(state%fitness_x)
            if (min_fitness < state%fitness_gbest) then
                j = minloc(state%fitness_x, dim=1)
                state%fitness_gbest = min_fitness
                state%x_gbest = state%x(:,j)
            end if
        end if

        relocate = .false.
        if (state%npart > 1) then
            select case (control%part_xchange)
            case (1)
                relocate = (state%futile_iter_count == maxval(state%futile_iter_count)) .and. &
                    (state%fitness_lbest >= maxval(state%fitness_lbest)) .and. &
                    (state%fitness_lbest < huge(1.0_dp))
            case (2)
                relocate = state%fitness_lbest > state%fitness_gbest
            case (3)
                idx = 0
                do i = 1, state%npart
                    if (state%fitness_lbest(i) <= state%fitness_gbest) cycle
                    if (idx == 0) then
                        idx = i
                    else if (state%futile_iter_count(i) > state%futile_iter_count(idx)) then
                        idx = i
                    end if
                end do
                if (idx > 0) relocate(idx) = .true.
            end select

            if (any(relocate)) then
                if (control%part_xchange == 1 .or. control%part_xchange == 2) then
                    do i = 1, state%npart
                        if (.not. relocate(i)) cycle
                        state%x_lbest(:,i) = state%x_gbest
                        state%fitness_lbest(i) = state%fitness_gbest
                        state%futile_iter_count(i) = 0
                    end do
                else if (control%part_xchange == 3) then
                    recent = minloc(state%futile_iter_count, dim=1)
                    do i = 1, state%npart
                        if (.not. relocate(i)) cycle
                        state%x_lbest(:,i) = state%x_lbest(:,recent)
                        state%fitness_lbest(i) = state%fitness_lbest(recent)
                        state%futile_iter_count(i) = state%futile_iter_count(recent)
                    end do
                end if
            end if
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

        if (state%init_function_calls+sum(state%function_calls) >= control%max_number_function_calls) &
            state%stop_code = PPSO_STOP_MAX_CALLS
        if (state%stop_code /= PPSO_STOP_NONE) return

        call propose_dds(control, state)
    end subroutine update_dds

    subroutine propose_dds(control, state)
        type(dds_control), intent(in) :: control
        type(dds_state), intent(inout) :: state
        real(dp), allocatable :: ranges(:)
        logical, allocatable :: selected(:)
        real(dp) :: p_inclusion
        integer :: i, j
        logical :: used_pending

        allocate(ranges(state%npar), selected(state%npar))
        ranges = control%r*(state%upper-state%lower)

        do i = 1, state%npart
            call append_pending(state%pending_initial, state%x(:,i), used_pending)
            if (used_pending) then
                state%v(:,i) = 0.0_dp
                cycle
            end if

            if (state%function_calls(i) <= 0) then
                p_inclusion = 1.0_dp
            else
                p_inclusion = 1.0_dp - log(real(state%function_calls(i),dp)) / &
                    log(max(real(state%main_max_calls,dp)/real(state%npart,dp), 1.0_dp+epsilon(1.0_dp)))
            end if
            selected = .false.
            do j = 1, state%npar
                selected(j) = state%rng%uniform() <= p_inclusion
            end do
            if (.not. any(selected)) selected(state%rng%randint(1,state%npar)) = .true.

            state%v(:,i) = 0.0_dp
            do j = 1, state%npar
                if (selected(j)) state%v(j,i) = state%rng%normal()*ranges(j)
            end do
            state%x(:,i) = state%x_lbest(:,i) + state%v(:,i)
            call reflect_bounds(state%x(:,i), state%lower, state%upper)
        end do
    end subroutine propose_dds

    subroutine reflect_bounds(x, lower, upper)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: lower(:), upper(:)
        logical, allocatable :: below(:), above(:), over(:)

        allocate(below(size(x)), above(size(x)), over(size(x)))
        below = x < lower
        above = x > upper
        where (below) x = lower + (lower-x)
        over = below .and. (x > upper)
        where (over) x = lower
        where (above) x = upper + (upper-x)
        over = above .and. (x < lower)
        where (over) x = upper
    end subroutine reflect_bounds

    subroutine sort_indices_by_value(values, order, first, last)
        real(dp), intent(in) :: values(:)
        integer, intent(inout) :: order(:)
        integer, intent(in) :: first, last
        integer :: i, j, key

        if (first > last) return
        do i = first+1, last
            key = order(i)
            j = i-1
            do while (j >= first)
                if (values(order(j)) <= values(key)) exit
                order(j+1) = order(j)
                j = j-1
            end do
            order(j+1) = key
        end do
    end subroutine sort_indices_by_value

end module ppso_dds
