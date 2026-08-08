! SPDX-License-Identifier: GPL-2.0-only
module soma_optimizer
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use soma_kinds, only : dp
    use soma_random, only : random_index, sample_without_replacement
    use soma_types, only : soma_bounds, soma_options, soma_result, soma_cost_function, &
                           strategy_all2one, strategy_t3a, strategy_pareto, all2one
    implicit none
    private
    public :: soma_optimize

contains

    subroutine evaluate_cost(callback, x, value, count)
        procedure(soma_cost_function) :: callback
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value
        integer, intent(inout) :: count

        value = callback(x)
        count = count + 1
    end subroutine evaluate_cost

    subroutine soma_optimize(cost_function, bnds, result, options, init)
        procedure(soma_cost_function) :: cost_function
        type(soma_bounds), intent(in) :: bnds
        type(soma_result), intent(out) :: result
        type(soma_options), intent(in), optional :: options
        real(dp), intent(in), optional :: init(:,:)

        type(soma_options) :: opt
        real(dp), allocatable :: population(:,:), costs(:), history_work(:)
        real(dp), allocatable :: steps(:), perturb(:,:), perturb_draw(:,:), direction(:), candidate(:)
        real(dp), allocatable :: random_steps(:), trial_costs(:,:)
        logical, allocatable :: to_migrate(:)
        integer, allocatable :: migrants(:), leader_pool(:), migrant_pool(:), order_idx(:)
        integer, allocatable :: eval_work(:), migrating_indices(:)
        integer :: n_params, pop_size, n_steps, migration_count, evaluation_count
        integer :: leader, n_migrants_selected, n_migrating, i, j, k, m, pidx, best_step
        integer :: leader_pool_count, migrant_pool_count, base_migrant_index
        integer :: random_flat_index
        real(dp) :: progress, perturbation_chance, leader_value, separation, sum_extremes
        real(dp) :: best_value, pi_value
        logical :: terminate_relative

        opt = all2one()
        if (present(options)) opt = options

        call validate_inputs(bnds, opt, init, result%status, result%message)
        if (result%status /= 0) return

        n_params = size(bnds%lower)
        pop_size = opt%population_size
        allocate(population(n_params,pop_size), costs(pop_size))
        allocate(history_work(0:opt%n_migrations), eval_work(0:opt%n_migrations))

        if (present(init)) then
            population = init
        else
            call random_number(population)
            do j = 1, pop_size
                population(:,j) = population(:,j) * (bnds%upper - bnds%lower) + bnds%lower
            end do
        end if

        evaluation_count = 0
        do j = 1, pop_size
            call evaluate_cost(cost_function, population(:,j), costs(j), evaluation_count)
        end do

        migration_count = 0
        history_work(0) = minval(costs)
        eval_work(0) = 0
        pi_value = acos(-1.0_dp)

        if (trim(opt%strategy) == strategy_all2one) then
            call make_all2one_steps(opt%path_length, opt%step_length, steps)
            n_steps = size(steps)
        else
            n_steps = opt%n_steps
            allocate(steps(n_steps))
        end if

        if (trim(opt%strategy) == strategy_pareto) then
            leader_pool_count = ceiling(real(pop_size,dp) / 25.0_dp)
            migrant_pool_count = ceiling(real(pop_size,dp) / 6.25_dp)
            base_migrant_index = int(floor(real(pop_size,dp) / 5.0_dp + 0.5_dp))
        else
            leader_pool_count = 0
            migrant_pool_count = 0
            base_migrant_index = 0
        end if

        do
            if (opt%n_migrations == 0) then
                progress = 0.0_dp
            else
                progress = real(migration_count,dp) / real(opt%n_migrations,dp)
            end if

            select case (trim(opt%strategy))
            case (strategy_all2one)
                leader = minloc(costs, dim=1)
                allocate(migrants(pop_size))
                migrants = [(i, i=1,pop_size)]
                n_migrants_selected = pop_size
                perturbation_chance = opt%perturbation_chance

            case (strategy_t3a)
                perturbation_chance = 0.05_dp + 0.9_dp * progress
                do k = 1, n_steps
                    steps(k) = real(k-1,dp) * (0.15_dp - 0.08_dp * progress)
                end do
                allocate(leader_pool(opt%leader_pool_size))
                call sample_without_replacement(pop_size, opt%leader_pool_size, leader_pool)
                leader = leader_pool(1)
                do i = 2, size(leader_pool)
                    if (costs(leader_pool(i)) < costs(leader)) leader = leader_pool(i)
                end do
                allocate(migrant_pool(opt%migrant_pool_size))
                call sample_without_replacement(pop_size, opt%migrant_pool_size, migrant_pool)
                call stable_order_subset(costs, migrant_pool, order_idx)
                n_migrants_selected = opt%n_migrants
                allocate(migrants(n_migrants_selected))
                migrants = order_idx(1:n_migrants_selected)
                deallocate(leader_pool, migrant_pool, order_idx)

            case (strategy_pareto)
                perturbation_chance = 0.5_dp + 0.45_dp * &
                    cos(opt%perturbation_frequency * pi_value * progress + pi_value)
                do k = 1, n_steps
                    steps(k) = real(k-1,dp) * (0.35_dp + 0.15_dp * &
                        cos(opt%step_frequency * pi_value * progress))
                end do
                call stable_order_all(costs, order_idx)
                leader = order_idx(random_index(leader_pool_count))
                n_migrants_selected = 1
                allocate(migrants(1))
                i = base_migrant_index + random_index(migrant_pool_count)
                migrants(1) = order_idx(i)
                deallocate(order_idx)

            case default
                result%status = 2
                result%message = "Unknown SOMA strategy"
                return
            end select

            leader_value = costs(leader)
            separation = maxval(costs) - minval(costs)
            sum_extremes = maxval(costs) + minval(costs)

            if (migration_count == opt%n_migrations) exit
            if (separation < opt%min_absolute_sep) exit
            terminate_relative = .false.
            if ((sum_extremes < 0.0_dp .or. sum_extremes > 0.0_dp) .and. &
                ieee_is_finite(separation) .and. ieee_is_finite(sum_extremes)) then
                terminate_relative = abs(separation / sum_extremes) < opt%min_relative_sep
            end if
            if (terminate_relative) exit

            allocate(perturb(n_params,pop_size), perturb_draw(n_params,n_migrants_selected))
            allocate(to_migrate(pop_size))
            perturb = 0.0_dp
            call random_number(perturb_draw)
            do m = 1, n_migrants_selected
                pidx = migrants(m)
                do i = 1, n_params
                    if (perturb_draw(i,m) < perturbation_chance) then
                        perturb(i,pidx) = 1.0_dp
                    else if (trim(opt%strategy) == strategy_pareto) then
                        perturb(i,pidx) = progress
                    else
                        perturb(i,pidx) = 0.0_dp
                    end if
                end do
            end do
            do pidx = 1, pop_size
                to_migrate(pidx) = sum(perturb(:,pidx)) > 0.0_dp
            end do
            to_migrate(leader) = .false.
            n_migrating = count(to_migrate)

            if (n_migrating == 0) then
                deallocate(migrants, perturb, perturb_draw, to_migrate)
                cycle
            end if

            allocate(migrating_indices(n_migrating))
            j = 0
            do pidx = 1, pop_size
                if (to_migrate(pidx)) then
                    j = j + 1
                    migrating_indices(j) = pidx
                end if
            end do

            ! The R implementation generates a full population_size x n_steps
            ! random repair array even though only n_migrating columns are used.
            ! Generate the same number of random values and use matching linear
            ! indices for out-of-bounds replacement.
            allocate(random_steps(n_params * pop_size * n_steps))
            call random_number(random_steps)
            do k = 1, n_steps
                do j = 1, n_migrating
                    m = migrating_indices(j)
                    do i = 1, n_params
                        random_flat_index = i + (j-1)*n_params + (k-1)*n_params*n_migrating
                        random_steps(random_flat_index) = random_steps(random_flat_index) * &
                            (bnds%upper(i)-bnds%lower(i)) + bnds%lower(i)
                    end do
                end do
            end do

            allocate(direction(n_params), candidate(n_params), trial_costs(n_migrating,n_steps))
            ! R's apply(populationSteps, 2:3, ...) traverses the migrant index
            ! fastest within each step. Evaluate all candidates before moving
            ! anyone, as the R implementation does.
            do k = 1, n_steps
                do j = 1, n_migrating
                    pidx = migrating_indices(j)
                    direction = population(:,pidx) - population(:,leader)
                    candidate = population(:,pidx) - steps(k) * direction * perturb(:,pidx)
                    do i = 1, n_params
                        if (candidate(i) < bnds%lower(i) .or. candidate(i) > bnds%upper(i)) then
                            random_flat_index = i + (j-1)*n_params + (k-1)*n_params*n_migrating
                            candidate(i) = random_steps(random_flat_index)
                        end if
                    end do
                    call evaluate_cost(cost_function, candidate, trial_costs(j,k), evaluation_count)
                end do
            end do

            do j = 1, n_migrating
                pidx = migrating_indices(j)
                direction = population(:,pidx) - population(:,leader)
                best_step = minloc(trial_costs(j,:), dim=1)
                best_value = trial_costs(j,best_step)
                candidate = population(:,pidx) - steps(best_step) * direction * perturb(:,pidx)
                do i = 1, n_params
                    if (candidate(i) < bnds%lower(i) .or. candidate(i) > bnds%upper(i)) then
                        ! Reconstruct the exact repair value selected for this step.
                        random_flat_index = i + (j-1)*n_params + &
                            (best_step-1)*n_params*n_migrating
                        candidate(i) = random_steps(random_flat_index)
                    end if
                end do
                population(:,pidx) = candidate
                costs(pidx) = best_value
            end do

            migration_count = migration_count + 1
            history_work(migration_count) = minval(costs)
            eval_work(migration_count) = evaluation_count

            deallocate(migrants, perturb, perturb_draw, to_migrate, migrating_indices)
            deallocate(random_steps, direction, candidate, trial_costs)
        end do

        leader = minloc(costs, dim=1)
        result%leader = leader
        result%population = population
        result%cost = costs
        result%history = history_work(0:migration_count)
        result%evaluations = eval_work(0:migration_count)
        result%migrations = migration_count
        result%status = 0
        result%message = "ok"

    end subroutine soma_optimize

    subroutine validate_inputs(bnds, opt, init, status, message)
        type(soma_bounds), intent(in) :: bnds
        type(soma_options), intent(in) :: opt
        real(dp), intent(in), optional :: init(:,:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: n

        status = 0
        message = "ok"
        if (.not. allocated(bnds%lower) .or. .not. allocated(bnds%upper)) then
            status = 1; message = "Bounds are not allocated"; return
        end if
        n = size(bnds%lower)
        if (n == 0 .or. size(bnds%upper) /= n) then
            status = 1; message = "Bounds must have equal nonzero length"; return
        end if
        if (any(.not. ieee_is_finite(bnds%lower)) .or. any(.not. ieee_is_finite(bnds%upper))) then
            status = 1; message = "Bounds must be finite"; return
        end if
        if (any(bnds%lower > bnds%upper)) then
            status = 1; message = "A lower bound exceeds its upper bound"; return
        end if
        if (opt%population_size < 2 .or. opt%n_migrations < 0) then
            status = 1; message = "Invalid population size or migration limit"; return
        end if
        select case (trim(opt%strategy))
        case (strategy_all2one)
            if (opt%step_length <= 0.0_dp .or. opt%path_length < 0.0_dp) then
                status = 1; message = "All2One requires nonnegative path length and positive step length"; return
            end if
        case (strategy_t3a)
            if (opt%n_steps < 1 .or. opt%leader_pool_size < 1 .or. opt%migrant_pool_size < 1 .or. &
                opt%leader_pool_size > opt%population_size .or. &
                opt%migrant_pool_size > opt%population_size .or. &
                opt%n_migrants < 1 .or. opt%n_migrants > opt%migrant_pool_size) then
                status = 1; message = "Invalid T3A pool or step sizes"; return
            end if
        case (strategy_pareto)
            if (opt%n_steps < 1) then
                status = 1; message = "Pareto requires at least one step"; return
            end if
        case default
            status = 1; message = "Unknown strategy"; return
        end select
        if (present(init)) then
            if (size(init,1) /= n .or. size(init,2) /= opt%population_size) then
                status = 1; message = "Initial population has the wrong shape"; return
            end if
        end if
    end subroutine validate_inputs

    subroutine make_all2one_steps(path_length, step_length, steps)
        real(dp), intent(in) :: path_length, step_length
        real(dp), allocatable, intent(out) :: steps(:)
        integer :: n, k
        real(dp) :: ratio

        ratio = path_length / step_length
        n = int(floor(ratio + 1.0e-10_dp)) + 1
        if (n < 1) n = 1
        allocate(steps(n))
        do k = 1, n
            steps(k) = real(k-1,dp) * step_length
        end do
    end subroutine make_all2one_steps

    subroutine stable_order_all(values, order_idx)
        real(dp), intent(in) :: values(:)
        integer, allocatable, intent(out) :: order_idx(:)
        integer :: i, j, key

        allocate(order_idx(size(values)))
        order_idx = [(i, i=1,size(values))]
        do i = 2, size(order_idx)
            key = order_idx(i)
            j = i - 1
            do while (j >= 1)
                if (values(order_idx(j)) <= values(key)) exit
                order_idx(j+1) = order_idx(j)
                j = j - 1
            end do
            order_idx(j+1) = key
        end do
    end subroutine stable_order_all

    subroutine stable_order_subset(values, subset, order_idx)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: subset(:)
        integer, allocatable, intent(out) :: order_idx(:)
        integer :: i, j, key

        allocate(order_idx(size(subset)))
        order_idx = subset
        do i = 2, size(order_idx)
            key = order_idx(i)
            j = i - 1
            do while (j >= 1)
                if (values(order_idx(j)) <= values(key)) exit
                order_idx(j+1) = order_idx(j)
                j = j - 1
            end do
            order_idx(j+1) = key
        end do
    end subroutine stable_order_subset

end module soma_optimizer
