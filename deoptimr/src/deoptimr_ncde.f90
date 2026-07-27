! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr_ncde
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use deoptimr_kinds, only: dp
    use deoptimr_interfaces, only: objective_function, constraint_function
    use deoptimr_rng, only: random_uniform, sample_without_replacement
    use deoptimr_types, only: ncde_control, ncde_result
    use deoptimr_utils, only: total_violation, transform_constraints, best_population_index
    use deoptimr_utils, only: handle_bounds, make_population, initialize_adaptation
    use deoptimr_utils, only: distance_squared, sort_indices_by_values, sort_indices_by_two_keys
    use deoptimr_utils, only: nearest_distances, almost_equal, nearest_neighbor_index
    implicit none
    private

    public :: ncde_optimize

contains

    subroutine ncde_optimize(lower, upper, objective, result, control, constraint, &
            n_constraints, n_equalities, equality_tolerance, initial_population)
        real(dp), intent(in) :: lower(:), upper(:)
        procedure(objective_function) :: objective
        type(ncde_result), intent(out) :: result
        type(ncde_control), intent(in), optional :: control
        procedure(constraint_function), optional :: constraint
        integer, intent(in), optional :: n_constraints, n_equalities
        real(dp), intent(in), optional :: equality_tolerance(:)
        real(dp), intent(in), optional :: initial_population(:, :)

        type(ncde_control) :: cfg
        real(dp), allocatable :: pop(:, :), pop_next(:, :), fpop(:), fpop_next(:)
        real(dp), allocatable :: f(:, :), f_next(:, :), cr(:), cr_next(:), pf(:), pf_next(:)
        real(dp), allocatable :: neighbors(:), neighbors_next(:)
        real(dp), allocatable :: hpop(:, :), hpop_next(:, :), violation(:), violation_next(:)
        real(dp), allocatable :: raw_h(:), htrial(:), trial(:), ftrial_vec(:)
        real(dp), allocatable :: distances(:), nearest(:), archive_x(:, :), archive_f(:)
        real(dp), allocatable :: archive_h(:, :)
        integer, allocatable :: distance_order(:), candidates(:), selected(:), population_order(:)
        integer, allocatable :: archive_order(:)
        integer :: d, np, np_random, fdim, ncon, meq, iteration
        integer :: i0, i, k, r(3), nneighbors, archive_count, best_index
        real(dp) :: ftrial, crtrial, pftrial, neighbor_trial, trial_violation
        real(dp) :: mu, ff, radius, best_fpop
        logical :: constrained, accepted, archive_candidate

        call prepare_ncde_control(lower, upper, control, cfg)
        d = size(lower)
        np_random = cfg%population_size
        call make_population(lower, upper, np_random, pop, initial_population)
        np = size(pop, 2)
        if (np < 4) error stop "ncde_optimize: population must contain at least four vectors"
        if (cfg%neighbor_lower == 0) cfg%neighbor_lower = max(3, np/20)
        if (cfg%neighbor_upper == 0) cfg%neighbor_upper = max(cfg%neighbor_lower, np/5)
        if (cfg%neighbor_lower < 3 .or. cfg%neighbor_upper > np - 1 .or. &
            cfg%neighbor_lower > cfg%neighbor_upper) then
            error stop "ncde_optimize: invalid neighborhood bounds"
        end if

        constrained = present(constraint)
        call constraint_dimensions(constrained, n_constraints, n_equalities, ncon, meq)
        fdim = merge(d, 1, cfg%use_jitter)

        allocate(pop_next(d, np), fpop(np), fpop_next(np), f(fdim, np), f_next(fdim, np))
        allocate(cr(np), cr_next(np), pf(np), pf_next(np), neighbors(np), neighbors_next(np))
        allocate(trial(d), ftrial_vec(d), distances(np - 1), distance_order(np - 1))
        allocate(candidates(np - 1), selected(3), nearest(np), population_order(np))
        call initialize_adaptation(d, np, cfg%f_lower, cfg%f_upper, cfg%cr_lower, cfg%cr_upper, &
            cfg%jitter_factor, cfg%use_jitter, f, cr, pf)
        do i = 1, np
            neighbors(i) = random_uniform(real(cfg%neighbor_lower, dp), real(cfg%neighbor_upper, dp))
            fpop(i) = objective(pop(:, i))
            if (ieee_is_nan(fpop(i))) error stop "ncde_optimize: objective returned NaN"
        end do

        if (constrained) then
            allocate(hpop(ncon, np), hpop_next(ncon, np), violation(np), violation_next(np))
            allocate(raw_h(ncon), htrial(ncon))
            do i = 1, np
                call constraint(pop(:, i), raw_h)
                if (any(ieee_is_nan(raw_h))) error stop "ncde_optimize: constraint returned NaN"
                call transform_constraints(raw_h, meq, equality_tolerance, hpop(:, i))
                violation(i) = total_violation(hpop(:, i))
            end do
            mu = median_value_local(violation)
            best_fpop = huge(1.0_dp)
        else
            mu = 0.0_dp
            best_fpop = minval(fpop)
        end if

        allocate(archive_x(d, max(cfg%archive_size, 1)), archive_f(max(cfg%archive_size, 1)))
        if (constrained) allocate(archive_h(ncon, max(cfg%archive_size, 1)))
        archive_count = 0
        if (cfg%niche_radius > 0.0_dp) then
            radius = cfg%niche_radius
        else
            radius = huge(1.0_dp)
        end if

        do iteration = 1, cfg%max_iterations
            if (cfg%niche_radius <= 0.0_dp) then
                call nearest_distances(pop, nearest)
                radius = min(radius, sum(nearest)/real(np, dp))
            end if

            pop_next = pop
            fpop_next = fpop
            f_next = f
            cr_next = cr
            pf_next = pf
            neighbors_next = neighbors
            if (constrained) then
                hpop_next = hpop
                violation_next = violation
            end if

            do i0 = 1, np
                i = mod(iteration + i0, np) + 1
                call draw_trial_adaptation(cfg, d, f(:, i), cr(i), pf(i), neighbors(i), &
                    ftrial_vec, crtrial, pftrial, neighbor_trial)

                call nearest_subpopulation(pop, i, distances, candidates, distance_order)
                nneighbors = min(np - 1, max(3, int(neighbor_trial)))
                call sample_without_replacement(candidates(distance_order(1:nneighbors)), 3, selected)
                r = selected
                call reproduce(pop(:, i), pop(:, r(1)), pop(:, r(2)), pop(:, r(3)), &
                    ftrial_vec, crtrial, pftrial, trial)
                call handle_bounds(trial, pop(:, r(1)), lower, upper)
                k = nearest_neighbor_index(trial, pop)

                accepted = .false.
                archive_candidate = .false.
                if (.not. constrained) then
                    ftrial = objective(trial)
                    if (ieee_is_nan(ftrial)) error stop "ncde_optimize: objective returned NaN"
                    if (ftrial <= fpop(k)) then
                        call set_trial_in_next(k)
                        accepted = .true.
                        archive_candidate = ftrial < best_fpop .or. &
                            almost_equal(ftrial, best_fpop, cfg%critical_tolerance)
                        if (ftrial < best_fpop) best_fpop = ftrial
                    end if
                else
                    call constraint(trial, raw_h)
                    if (any(ieee_is_nan(raw_h))) error stop "ncde_optimize: constraint returned NaN"
                    call transform_constraints(raw_h, meq, equality_tolerance, htrial)
                    trial_violation = total_violation(htrial)
                    if (trial_violation > mu) then
                        if (trial_violation <= violation(k)) then
                            call set_constrained_trial_in_next(k, .false.)
                            accepted = .true.
                        end if
                    else if (violation(k) > mu) then
                        ftrial = objective(trial)
                        if (ieee_is_nan(ftrial)) error stop "ncde_optimize: objective returned NaN"
                        call set_constrained_trial_in_next(k, .true.)
                        accepted = .true.
                        if (meq == 0) then
                            ff = real(count(violation <= mu), dp)/real(np, dp)
                            mu = mu*(1.0_dp - ff/real(np, dp))
                        end if
                    else
                        ftrial = objective(trial)
                        if (ieee_is_nan(ftrial)) error stop "ncde_optimize: objective returned NaN"
                        if (ftrial <= fpop(k)) then
                            call set_constrained_trial_in_next(k, .true.)
                            accepted = .true.
                            ff = real(count(violation <= mu), dp)/real(np, dp)
                            mu = mu*(1.0_dp - ff/real(np, dp))
                        end if
                    end if
                    if (accepted .and. all(htrial <= 0.0_dp)) then
                        archive_candidate = ftrial < best_fpop .or. &
                            almost_equal(ftrial, best_fpop, cfg%critical_tolerance)
                        if (ftrial < best_fpop) best_fpop = ftrial
                    end if
                end if

                if (archive_candidate) then
                    call update_archive(trial, ftrial, htrial, constrained, radius, cfg%archive_size, &
                        archive_x, archive_f, archive_h, archive_count)
                    if (cfg%reinitialize_archive_neighbors .and. archive_count > 0) then
                        call reinitialize_member(k)
                    end if
                end if
            end do

            pop = pop_next
            fpop = fpop_next
            f = f_next
            cr = cr_next
            pf = pf_next
            neighbors = neighbors_next
            if (constrained) then
                hpop = hpop_next
                violation = violation_next
            end if

            if (cfg%trace .and. mod(iteration, cfg%trace_interval) == 0) then
                if (constrained) then
                    best_index = best_population_index(fpop, violation, mu, .true.)
                else
                    best_index = best_population_index(fpop, constrained=.false.)
                end if
                write(*, '(i0,2x,a,es13.5,2x,a,i0,2x,a,es13.5)') iteration, &
                    "radius=", radius, "archive=", archive_count, "best=", fpop(best_index)
            end if
        end do

        result%iterations = cfg%max_iterations
        result%final_niche_radius = radius
        if (archive_count > 0) then
            allocate(archive_order(archive_count))
            call sort_indices_by_values(archive_f(1:archive_count), archive_order)
            allocate(result%solution_archive(d, archive_count), result%objective_archive(archive_count))
            result%solution_archive = archive_x(:, archive_order)
            result%objective_archive = archive_f(archive_order)
            if (constrained) then
                allocate(result%constraint_archive(ncon, archive_count))
                result%constraint_archive = archive_h(:, archive_order)
            end if
        else
            allocate(result%solution_archive(d, 0), result%objective_archive(0))
            if (constrained) allocate(result%constraint_archive(ncon, 0))
        end if

        if (constrained) then
            call sort_population_constrained(hpop, fpop, population_order)
        else
            call sort_indices_by_values(fpop, population_order)
        end if
        allocate(result%solution_population(d, np), result%objective_population(np))
        result%solution_population = pop(:, population_order)
        result%objective_population = fpop(population_order)
        if (constrained) then
            allocate(result%constraint_population(ncon, np), result%total_violation_population(np))
            result%constraint_population = hpop(:, population_order)
            result%total_violation_population = violation(population_order)
        end if

    contains

        subroutine set_trial_in_next(index)
            integer, intent(in) :: index
            pop_next(:, index) = trial
            fpop_next(index) = ftrial
            if (cfg%use_jitter) then
                f_next(:, index) = ftrial_vec
            else
                f_next(1, index) = ftrial_vec(1)
            end if
            cr_next(index) = crtrial
            pf_next(index) = pftrial
            neighbors_next(index) = neighbor_trial
        end subroutine set_trial_in_next

        subroutine set_constrained_trial_in_next(index, set_objective)
            integer, intent(in) :: index
            logical, intent(in) :: set_objective
            call set_trial_in_next(index)
            if (.not. set_objective) fpop_next(index) = fpop(index)
            hpop_next(:, index) = htrial
            violation_next(index) = trial_violation
        end subroutine set_constrained_trial_in_next

        subroutine reinitialize_member(index)
            integer, intent(in) :: index
            integer :: m
            do m = 1, d
                pop_next(m, index) = random_uniform(lower(m), upper(m))
            end do
            fpop_next(index) = objective(pop_next(:, index))
            if (ieee_is_nan(fpop_next(index))) error stop "ncde_optimize: objective returned NaN"
            call draw_f(cfg, d, ftrial_vec)
            if (cfg%use_jitter) then
                f_next(:, index) = ftrial_vec
            else
                f_next(1, index) = ftrial_vec(1)
            end if
            cr_next(index) = random_uniform(cfg%cr_lower, cfg%cr_upper)
            pf_next(index) = random_uniform()
            neighbors_next(index) = random_uniform(real(cfg%neighbor_lower, dp), &
                real(cfg%neighbor_upper, dp))
            if (constrained) then
                call constraint(pop_next(:, index), raw_h)
                call transform_constraints(raw_h, meq, equality_tolerance, hpop_next(:, index))
                violation_next(index) = total_violation(hpop_next(:, index))
            end if
        end subroutine reinitialize_member
    end subroutine ncde_optimize

    subroutine prepare_ncde_control(lower, upper, control, cfg)
        real(dp), intent(in) :: lower(:), upper(:)
        type(ncde_control), intent(in), optional :: control
        type(ncde_control), intent(out) :: cfg

        if (size(lower) <= 0 .or. size(upper) /= size(lower)) error stop "NCDE: invalid bounds"
        if (any(lower > upper)) error stop "NCDE: lower bound exceeds upper"
        cfg = ncde_control()
        if (present(control)) cfg = control
        if (cfg%population_size < 0 .or. cfg%archive_size < 0 .or. cfg%max_iterations < 0) then
            error stop "NCDE: negative size"
        end if
        if (cfg%critical_tolerance <= 0.0_dp) error stop "NCDE: critical tolerance must be positive"
        if (cfg%f_lower > cfg%f_upper .or. cfg%cr_lower < 0.0_dp .or. &
            cfg%cr_lower > 1.0_dp .or. cfg%cr_upper < cfg%cr_lower) then
            error stop "NCDE: invalid adaptation interval"
        end if
        if (cfg%tau_f < 0.0_dp .or. cfg%tau_f > 1.0_dp .or. &
            cfg%tau_cr < 0.0_dp .or. cfg%tau_cr > 1.0_dp .or. &
            cfg%tau_pf < 0.0_dp .or. cfg%tau_pf > 1.0_dp .or. &
            cfg%tau_neighbors < 0.0_dp .or. cfg%tau_neighbors > 1.0_dp) then
            error stop "NCDE: adaptation probabilities must be in [0,1]"
        end if
        if (cfg%trace_interval < 1) error stop "NCDE: trace interval must be positive"
    end subroutine prepare_ncde_control

    subroutine constraint_dimensions(constrained, n_constraints, n_equalities, ncon, meq)
        logical, intent(in) :: constrained
        integer, intent(in), optional :: n_constraints, n_equalities
        integer, intent(out) :: ncon, meq

        if (constrained) then
            if (.not. present(n_constraints)) error stop "NCDE: n_constraints is required"
            ncon = n_constraints
            meq = 0
            if (present(n_equalities)) meq = n_equalities
            if (ncon <= 0 .or. meq < 0 .or. meq > ncon) error stop "NCDE: invalid constraints"
        else
            ncon = 0
            meq = 0
        end if
    end subroutine constraint_dimensions

    subroutine draw_trial_adaptation(cfg, d, f_current, cr_current, pf_current, neighbor_current, &
            ftrial, crtrial, pftrial, neighbor_trial)
        type(ncde_control), intent(in) :: cfg
        integer, intent(in) :: d
        real(dp), intent(in) :: f_current(:), cr_current, pf_current, neighbor_current
        real(dp), intent(out) :: ftrial(d), crtrial, pftrial, neighbor_trial

        if (random_uniform() <= cfg%tau_f) then
            call draw_f(cfg, d, ftrial)
        else if (cfg%use_jitter) then
            ftrial = f_current
        else
            ftrial = f_current(1)
        end if
        if (random_uniform() <= cfg%tau_cr) then
            crtrial = random_uniform(cfg%cr_lower, cfg%cr_upper)
        else
            crtrial = cr_current
        end if
        if (random_uniform() <= cfg%tau_pf) then
            pftrial = random_uniform()
        else
            pftrial = pf_current
        end if
        if (random_uniform() <= cfg%tau_neighbors) then
            neighbor_trial = random_uniform(real(cfg%neighbor_lower, dp), real(cfg%neighbor_upper, dp))
        else
            neighbor_trial = neighbor_current
        end if
    end subroutine draw_trial_adaptation

    subroutine draw_f(cfg, d, ftrial)
        type(ncde_control), intent(in) :: cfg
        integer, intent(in) :: d
        real(dp), intent(out) :: ftrial(d)
        real(dp) :: scale
        integer :: j

        scale = random_uniform(cfg%f_lower, cfg%f_upper)
        if (cfg%use_jitter) then
            do j = 1, d
                ftrial(j) = scale*(1.0_dp + cfg%jitter_factor*random_uniform(-0.5_dp, 0.5_dp))
            end do
        else
            ftrial = scale
        end if
    end subroutine draw_f

    subroutine reproduce(target, base, r1, r2, ftrial, crtrial, pftrial, trial)
        use deoptimr_rng, only: random_integer
        real(dp), intent(in) :: target(:), base(:), r1(:), r2(:), ftrial(:)
        real(dp), intent(in) :: crtrial, pftrial
        real(dp), intent(out) :: trial(size(target))
        logical, allocatable :: ignore(:)
        integer :: j

        allocate(ignore(size(target)))
        do j = 1, size(target)
            ignore(j) = random_uniform() > crtrial
        end do
        if (all(ignore)) ignore(random_integer(1, size(target))) = .false.
        if (random_uniform() <= pftrial) then
            trial = base + ftrial*(r1 - r2)
        else
            trial = base + 0.5_dp*(ftrial + 1.0_dp)*(r1 + r2 - 2.0_dp*base)
        end if
        where (ignore) trial = target
    end subroutine reproduce

    subroutine nearest_subpopulation(population, excluded, distances, candidates, order)
        real(dp), intent(in) :: population(:, :)
        integer, intent(in) :: excluded
        real(dp), intent(out) :: distances(:)
        integer, intent(out) :: candidates(:), order(:)
        integer :: j, k

        k = 0
        do j = 1, size(population, 2)
            if (j == excluded) cycle
            k = k + 1
            candidates(k) = j
            distances(k) = sqrt(distance_squared(population(:, excluded), population(:, j)))
        end do
        call sort_indices_by_values(distances, order)
    end subroutine nearest_subpopulation

    subroutine update_archive(x, value, h, constrained, radius, capacity, &
            archive_x, archive_f, archive_h, count_archive)
        real(dp), intent(in) :: x(:), value
        real(dp), intent(in), optional :: h(:)
        logical, intent(in) :: constrained
        real(dp), intent(in) :: radius
        integer, intent(in) :: capacity
        real(dp), intent(inout) :: archive_x(:, :), archive_f(:)
        real(dp), intent(inout), optional :: archive_h(:, :)
        integer, intent(inout) :: count_archive
        integer :: j, nearest_index
        real(dp) :: current, best_distance

        if (capacity == 0) return
        if (count_archive == 0) then
            count_archive = 1
            archive_x(:, 1) = x
            archive_f(1) = value
            if (constrained) archive_h(:, 1) = h
            return
        end if

        nearest_index = 0
        best_distance = huge(1.0_dp)
        do j = 1, count_archive
            current = sqrt(distance_squared(x, archive_x(:, j)))
            if (current <= radius .and. current < best_distance) then
                best_distance = current
                nearest_index = j
            end if
        end do
        if (nearest_index > 0) then
            if (value < archive_f(nearest_index)) then
                archive_x(:, nearest_index) = x
                archive_f(nearest_index) = value
                if (constrained) archive_h(:, nearest_index) = h
            end if
        else if (count_archive < capacity) then
            count_archive = count_archive + 1
            archive_x(:, count_archive) = x
            archive_f(count_archive) = value
            if (constrained) archive_h(:, count_archive) = h
        end if
    end subroutine update_archive

    subroutine sort_population_constrained(h, cost, order)
        real(dp), intent(in) :: h(:, :), cost(:)
        integer, intent(out) :: order(size(cost))
        real(dp), allocatable :: infeasible(:)
        integer :: j

        allocate(infeasible(size(cost)))
        do j = 1, size(cost)
            infeasible(j) = merge(1.0_dp, 0.0_dp, any(h(:, j) > 0.0_dp))
        end do
        call sort_indices_by_two_keys(infeasible, cost, order)
    end subroutine sort_population_constrained

    function median_value_local(x) result(value)
        use deoptimr_utils, only: median_value
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = median_value(x)
    end function median_value_local
end module deoptimr_ncde
