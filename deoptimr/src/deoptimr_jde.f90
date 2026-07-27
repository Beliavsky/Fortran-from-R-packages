! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr_jde
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use deoptimr_kinds, only: dp
    use deoptimr_interfaces, only: objective_function, constraint_function
    use deoptimr_rng, only: random_uniform, random_integer
    use deoptimr_types, only: jde_control, de_result
    use deoptimr_utils, only: objective_spread, total_violation, transform_constraints
    use deoptimr_utils, only: best_population_index, handle_bounds, make_population
    use deoptimr_utils, only: initialize_adaptation, choose_distinct_indices
    implicit none
    private

    public :: jde_optimize, spjde_optimize

contains

    subroutine jde_optimize(lower, upper, objective, result, control, constraint, &
            n_constraints, n_equalities, equality_tolerance, initial_population)
        real(dp), intent(in) :: lower(:), upper(:)
        procedure(objective_function) :: objective
        type(de_result), intent(out) :: result
        type(jde_control), intent(in), optional :: control
        procedure(constraint_function), optional :: constraint
        integer, intent(in), optional :: n_constraints, n_equalities
        real(dp), intent(in), optional :: equality_tolerance(:)
        real(dp), intent(in), optional :: initial_population(:, :)

        type(jde_control) :: cfg
        real(dp), allocatable :: pop(:, :), fpop(:), f(:, :), cr(:), pf(:)
        real(dp), allocatable :: hpop(:, :), violation(:), raw_h(:), htrial(:)
        real(dp), allocatable :: trial(:), ftrial_vec(:)
        integer :: d, np_random, np, fdim, maxiter, iteration, convergence
        integer :: i0, i, r(3), best, ncon, meq
        real(dp) :: ftrial, crtrial, pftrial, mu, convergence_measure
        logical :: constrained, accepted, continue_search

        call prepare_jde_control(lower, upper, control, cfg, .false.)
        d = size(lower)
        np_random = cfg%population_size
        call make_population(lower, upper, np_random, pop, initial_population)
        np = size(pop, 2)
        if (np < 4) error stop "jde_optimize: population must contain at least four vectors"
        maxiter = cfg%max_iterations

        constrained = present(constraint)
        call constraint_dimensions(constrained, n_constraints, n_equalities, ncon, meq)

        fdim = merge(d, 1, cfg%use_jitter)
        allocate(f(fdim, np), cr(np), pf(np), fpop(np), trial(d), ftrial_vec(d))
        call initialize_adaptation(d, np, cfg%f_lower, cfg%f_upper, 0.0_dp, 1.0_dp, &
            cfg%jitter_factor, cfg%use_jitter, f, cr, pf)
        do i = 1, np
            fpop(i) = objective(pop(:, i))
            if (ieee_is_nan(fpop(i))) error stop "jde_optimize: objective returned NaN"
        end do

        if (constrained) then
            allocate(hpop(ncon, np), violation(np), raw_h(ncon), htrial(ncon))
            do i = 1, np
                call constraint(pop(:, i), raw_h)
                if (any(ieee_is_nan(raw_h))) error stop "jde_optimize: constraint returned NaN"
                call transform_constraints(raw_h, meq, equality_tolerance, hpop(:, i))
                violation(i) = total_violation(hpop(:, i))
            end do
            mu = median_local(violation)
            best = best_population_index(fpop, violation, mu, .true.)
        else
            mu = 0.0_dp
            best = best_population_index(fpop, constrained=.false.)
        end if

        convergence_measure = objective_spread(fpop, fpop(best), cfg%compare_to, cfg%objective_scale)
        iteration = 0
        convergence = 0
        if (constrained) then
            continue_search = continuation_rule(convergence_measure, cfg%tolerance, .true., hpop, best)
        else
            continue_search = continuation_rule(convergence_measure, cfg%tolerance, .false., best=best)
        end if

        do while (continue_search)
            if (iteration >= maxiter) then
                convergence = 1
                exit
            end if
            iteration = iteration + 1

            do i0 = 1, np
                i = mod(iteration + i0, np) + 1
                call trial_adaptation(cfg, d, f(:, i), cr(i), pf(i), ftrial_vec, crtrial, pftrial)
                call choose_distinct_indices(np, i, r)
                call reproduce(pop(:, i), pop(:, r(1)), pop(:, r(2)), pop(:, r(3)), &
                    ftrial_vec, crtrial, pftrial, trial)
                call handle_bounds(trial, pop(:, r(1)), lower, upper)

                accepted = .false.
                if (.not. constrained) then
                    ftrial = objective(trial)
                    if (ieee_is_nan(ftrial)) error stop "jde_optimize: objective returned NaN"
                    if (ftrial <= fpop(i)) accepted = .true.
                else
                    call constraint(trial, raw_h)
                    if (any(ieee_is_nan(raw_h))) error stop "jde_optimize: constraint returned NaN"
                    call transform_constraints(raw_h, meq, equality_tolerance, htrial)
                    call constrained_selection_jde(i, trial, htrial, objective, pop, fpop, hpop, &
                        violation, mu, meq, accepted, ftrial)
                end if

                if (accepted) then
                    if (.not. constrained) then
                        pop(:, i) = trial
                        fpop(i) = ftrial
                    end if
                    if (cfg%use_jitter) then
                        f(:, i) = ftrial_vec
                    else
                        f(1, i) = ftrial_vec(1)
                    end if
                    cr(i) = crtrial
                    pf(i) = pftrial
                end if

                if (constrained) then
                    best = best_population_index(fpop, violation, mu, .true.)
                else
                    best = best_population_index(fpop, constrained=.false.)
                end if
            end do

            convergence_measure = objective_spread(fpop, fpop(best), cfg%compare_to, cfg%objective_scale)
            if (cfg%trace .and. mod(iteration, cfg%trace_interval) == 0) then
                write(*, '(i0,2x,a,es13.5,2x,a,es13.5)') iteration, "spread=", &
                    convergence_measure, "best=", fpop(best)
            end if
            if (constrained) then
                continue_search = continuation_rule(convergence_measure, cfg%tolerance, .true., hpop, best)
            else
                continue_search = continuation_rule(convergence_measure, cfg%tolerance, .false., best=best)
            end if
        end do

        allocate(result%parameters(d))
        result%parameters = pop(:, best)
        result%value = fpop(best)
        result%iterations = iteration
        result%convergence = convergence
        if (cfg%save_population) then
            allocate(result%population(d, np), result%population_cost(np))
            result%population = pop
            result%population_cost = fpop
            if (constrained) then
                allocate(result%population_constraints(ncon, np), result%total_violation(np))
                result%population_constraints = hpop
                result%total_violation = violation
            end if
        end if
    end subroutine jde_optimize

    subroutine spjde_optimize(lower, upper, objective, result, control, constraint, &
            n_constraints, n_equalities, equality_tolerance, initial_population)
        real(dp), intent(in) :: lower(:), upper(:)
        procedure(objective_function) :: objective
        type(de_result), intent(out) :: result
        type(jde_control), intent(in), optional :: control
        procedure(constraint_function), optional :: constraint
        integer, intent(in), optional :: n_constraints, n_equalities
        real(dp), intent(in), optional :: equality_tolerance(:)
        real(dp), intent(in), optional :: initial_population(:, :)

        type(jde_control) :: cfg
        real(dp), allocatable :: pop(:, :), trial(:, :), fpop(:), ftrial(:)
        real(dp), allocatable :: f(:, :), ftrial_adapt(:, :), cr(:), crtrial(:), pf(:), pftrial(:)
        real(dp), allocatable :: hpop(:, :), htrial(:, :), violation(:), violation_trial(:), raw_h(:)
        real(dp), allocatable :: fvec(:)
        integer :: d, np, np_random, fdim, maxiter, iteration, convergence
        integer :: i0, i, r(3), best, ncon, meq, naccepted
        real(dp) :: mu, ff, convergence_measure
        logical :: constrained, continue_search

        call prepare_jde_control(lower, upper, control, cfg, .true.)
        d = size(lower)
        np_random = cfg%population_size
        call make_population(lower, upper, np_random, pop, initial_population)
        np = size(pop, 2)
        if (np < 4) error stop "spjde_optimize: population must contain at least four vectors"
        maxiter = cfg%max_iterations
        constrained = present(constraint)
        call constraint_dimensions(constrained, n_constraints, n_equalities, ncon, meq)

        fdim = merge(d, 1, cfg%use_jitter)
        allocate(trial(d, np), fpop(np), ftrial(np), f(fdim, np), ftrial_adapt(fdim, np))
        allocate(cr(np), crtrial(np), pf(np), pftrial(np), fvec(d))
        call initialize_adaptation(d, np, cfg%f_lower, cfg%f_upper, cfg%cr_lower, cfg%cr_upper, &
            cfg%jitter_factor, cfg%use_jitter, f, cr, pf)
        ftrial_adapt = f
        crtrial = cr
        pftrial = pf
        trial = pop

        if (constrained) then
            allocate(hpop(ncon, np), htrial(ncon, np), violation(np), violation_trial(np), raw_h(ncon))
            do i = 1, np
                call constraint(pop(:, i), raw_h)
                if (any(ieee_is_nan(raw_h))) error stop "spjde_optimize: constraint returned NaN"
                call transform_constraints(raw_h, meq, equality_tolerance, hpop(:, i))
                violation(i) = total_violation(hpop(:, i))
            end do
            mu = median_local(violation)
            fpop = huge(1.0_dp)
            do i = 1, np
                if (violation(i) <= mu) fpop(i) = objective(pop(:, i))
                if (ieee_is_nan(fpop(i))) error stop "spjde_optimize: objective returned NaN"
            end do
            best = best_population_index(fpop, violation, mu, .true.)
        else
            do i = 1, np
                fpop(i) = objective(pop(:, i))
                if (ieee_is_nan(fpop(i))) error stop "spjde_optimize: objective returned NaN"
            end do
            mu = 0.0_dp
            best = best_population_index(fpop, constrained=.false.)
        end if

        convergence_measure = objective_spread(fpop, fpop(best), cfg%compare_to, cfg%objective_scale)
        iteration = 0
        convergence = 0
        if (constrained) then
            continue_search = continuation_rule(convergence_measure, cfg%tolerance, .true., hpop, best)
        else
            continue_search = continuation_rule(convergence_measure, cfg%tolerance, .false., best=best)
        end if

        do while (continue_search)
            if (iteration >= maxiter) then
                convergence = 1
                exit
            end if
            iteration = iteration + 1

            do i = 1, np
                if (random_uniform() <= cfg%tau_f) then
                    call draw_f(cfg, d, fvec)
                    if (cfg%use_jitter) then
                        ftrial_adapt(:, i) = fvec
                    else
                        ftrial_adapt(1, i) = fvec(1)
                    end if
                end if
                if (random_uniform() <= cfg%tau_cr) crtrial(i) = random_uniform(cfg%cr_lower, cfg%cr_upper)
                if (random_uniform() <= cfg%tau_pf) pftrial(i) = random_uniform()
            end do

            do i0 = 1, np
                i = mod(iteration + i0, np) + 1
                call choose_distinct_indices(np, i, r)
                if (cfg%use_jitter) then
                    fvec = ftrial_adapt(:, i)
                else
                    fvec = ftrial_adapt(1, i)
                end if
                call reproduce(pop(:, i), pop(:, r(1)), pop(:, r(2)), pop(:, r(3)), &
                    fvec, crtrial(i), pftrial(i), trial(:, i))
                call handle_bounds(trial(:, i), pop(:, r(1)), lower, upper)
            end do

            ftrial = huge(1.0_dp)
            if (constrained) then
                do i = 1, np
                    call constraint(trial(:, i), raw_h)
                    if (any(ieee_is_nan(raw_h))) error stop "spjde_optimize: constraint returned NaN"
                    call transform_constraints(raw_h, meq, equality_tolerance, htrial(:, i))
                    violation_trial(i) = total_violation(htrial(:, i))
                    if (violation_trial(i) <= mu) then
                        ftrial(i) = objective(trial(:, i))
                        if (ieee_is_nan(ftrial(i))) error stop "spjde_optimize: objective returned NaN"
                    end if
                end do
                naccepted = 0
                do i = 1, np
                    if (violation_trial(i) > mu) then
                        if (violation_trial(i) <= violation(i)) then
                            call accept_spjde(i, .false.)
                        end if
                    else if (violation(i) > mu) then
                        call accept_spjde(i, .true.)
                        if (meq == 0) naccepted = naccepted + 1
                    else if (ftrial(i) <= fpop(i)) then
                        call accept_spjde(i, .true.)
                        naccepted = naccepted + 1
                    end if
                end do
                if (naccepted > 0) then
                    ff = real(count(violation <= mu), dp)/real(np, dp)
                    mu = mu*(1.0_dp - ff/real(np, dp))**naccepted
                end if
                best = best_population_index(fpop, violation, mu, .true.)
            else
                do i = 1, np
                    ftrial(i) = objective(trial(:, i))
                    if (ieee_is_nan(ftrial(i))) error stop "spjde_optimize: objective returned NaN"
                    if (ftrial(i) <= fpop(i)) call accept_spjde(i, .true.)
                end do
                best = best_population_index(fpop, constrained=.false.)
            end if

            convergence_measure = objective_spread(fpop, fpop(best), cfg%compare_to, cfg%objective_scale)
            if (cfg%trace .and. mod(iteration, cfg%trace_interval) == 0) then
                write(*, '(i0,2x,a,es13.5,2x,a,es13.5)') iteration, "spread=", &
                    convergence_measure, "best=", fpop(best)
            end if
            if (constrained) then
                continue_search = continuation_rule(convergence_measure, cfg%tolerance, .true., hpop, best)
            else
                continue_search = continuation_rule(convergence_measure, cfg%tolerance, .false., best=best)
            end if
        end do

        allocate(result%parameters(d))
        result%parameters = pop(:, best)
        result%value = fpop(best)
        result%iterations = iteration
        result%convergence = convergence
        if (cfg%save_population) then
            allocate(result%population(d, np), result%population_cost(np))
            result%population = pop
            result%population_cost = fpop
            if (constrained) then
                allocate(result%population_constraints(ncon, np), result%total_violation(np))
                result%population_constraints = hpop
                result%total_violation = violation
            end if
        end if

    contains

        subroutine accept_spjde(index, evaluate_objective)
            integer, intent(in) :: index
            logical, intent(in) :: evaluate_objective

            pop(:, index) = trial(:, index)
            if (evaluate_objective) fpop(index) = ftrial(index)
            if (constrained) then
                hpop(:, index) = htrial(:, index)
                violation(index) = violation_trial(index)
            end if
            f(:, index) = ftrial_adapt(:, index)
            cr(index) = crtrial(index)
            pf(index) = pftrial(index)
        end subroutine accept_spjde
    end subroutine spjde_optimize

    subroutine prepare_jde_control(lower, upper, control, cfg, synchronous)
        real(dp), intent(in) :: lower(:), upper(:)
        type(jde_control), intent(in), optional :: control
        type(jde_control), intent(out) :: cfg
        logical, intent(in) :: synchronous
        integer :: d

        d = size(lower)
        if (d <= 0 .or. size(upper) /= d) error stop "DE optimizer: invalid bounds"
        if (any(lower > upper)) error stop "DE optimizer: lower bound exceeds upper"
        cfg = jde_control()
        if (present(control)) cfg = control
        if (cfg%population_size == 0) cfg%population_size = 10*d
        if (cfg%max_iterations == 0) cfg%max_iterations = 2000*d
        if (cfg%population_size < 0 .or. cfg%max_iterations < 0) error stop "DE optimizer: negative size"
        if (cfg%f_lower > cfg%f_upper) error stop "DE optimizer: invalid F interval"
        if (cfg%cr_lower < 0.0_dp .or. cfg%cr_lower > 1.0_dp .or. &
            cfg%cr_upper < cfg%cr_lower) error stop "DE optimizer: invalid CR interval"
        if (.not. synchronous) then
            cfg%cr_lower = 0.0_dp
            cfg%cr_upper = 1.0_dp
        end if
        if (cfg%tau_f < 0.0_dp .or. cfg%tau_f > 1.0_dp .or. &
            cfg%tau_cr < 0.0_dp .or. cfg%tau_cr > 1.0_dp .or. &
            cfg%tau_pf < 0.0_dp .or. cfg%tau_pf > 1.0_dp) then
            error stop "DE optimizer: adaptation probabilities must be in [0,1]"
        end if
        if (cfg%trace_interval < 1) error stop "DE optimizer: trace interval must be positive"
        if (cfg%objective_scale <= 0.0_dp) error stop "DE optimizer: objective scale must be positive"
        if (trim(cfg%compare_to) /= "median" .and. trim(cfg%compare_to) /= "max") then
            error stop "DE optimizer: compare_to must be median or max"
        end if
    end subroutine prepare_jde_control

    subroutine constraint_dimensions(constrained, n_constraints, n_equalities, ncon, meq)
        logical, intent(in) :: constrained
        integer, intent(in), optional :: n_constraints, n_equalities
        integer, intent(out) :: ncon, meq

        if (constrained) then
            if (.not. present(n_constraints)) error stop "DE optimizer: n_constraints is required"
            ncon = n_constraints
            meq = 0
            if (present(n_equalities)) meq = n_equalities
            if (ncon <= 0 .or. meq < 0 .or. meq > ncon) error stop "DE optimizer: invalid constraints"
        else
            ncon = 0
            meq = 0
            if (present(n_constraints)) then
                if (n_constraints /= 0) error stop "DE optimizer: constraint callback is absent"
            end if
        end if
    end subroutine constraint_dimensions

    subroutine trial_adaptation(cfg, d, f_current, cr_current, pf_current, ftrial, crtrial, pftrial)
        type(jde_control), intent(in) :: cfg
        integer, intent(in) :: d
        real(dp), intent(in) :: f_current(:), cr_current, pf_current
        real(dp), intent(out) :: ftrial(d), crtrial, pftrial

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
    end subroutine trial_adaptation

    subroutine draw_f(cfg, d, ftrial)
        type(jde_control), intent(in) :: cfg
        integer, intent(in) :: d
        real(dp), intent(out) :: ftrial(d)
        integer :: j
        real(dp) :: scale

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
        real(dp), intent(in) :: target(:), base(:), r1(:), r2(:), ftrial(:)
        real(dp), intent(in) :: crtrial, pftrial
        real(dp), intent(out) :: trial(size(target))
        logical, allocatable :: ignore(:)
        integer :: j, forced

        allocate(ignore(size(target)))
        do j = 1, size(target)
            ignore(j) = random_uniform() > crtrial
        end do
        if (all(ignore)) then
            forced = random_integer(1, size(target))
            ignore(forced) = .false.
        end if
        if (random_uniform() <= pftrial) then
            trial = base + ftrial*(r1 - r2)
        else
            trial = base + 0.5_dp*(ftrial + 1.0_dp)*(r1 + r2 - 2.0_dp*base)
        end if
        where (ignore) trial = target
    end subroutine reproduce

    subroutine constrained_selection_jde(i, trial, htrial, objective, pop, fpop, hpop, &
            violation, mu, meq, accepted, ftrial)
        integer, intent(in) :: i
        real(dp), intent(in) :: trial(:), htrial(:)
        procedure(objective_function) :: objective
        real(dp), intent(inout) :: pop(:, :), fpop(:), hpop(:, :), violation(:), mu
        integer, intent(in) :: meq
        logical, intent(out) :: accepted
        real(dp), intent(out) :: ftrial
        real(dp) :: trial_violation, ff

        accepted = .false.
        ftrial = huge(1.0_dp)
        trial_violation = total_violation(htrial)
        if (trial_violation > mu) then
            if (trial_violation <= violation(i)) then
                pop(:, i) = trial
                hpop(:, i) = htrial
                violation(i) = trial_violation
                accepted = .true.
            end if
        else if (violation(i) > mu) then
            ftrial = objective(trial)
            if (ieee_is_nan(ftrial)) error stop "jde_optimize: objective returned NaN"
            pop(:, i) = trial
            fpop(i) = ftrial
            hpop(:, i) = htrial
            violation(i) = trial_violation
            accepted = .true.
            if (meq == 0) then
                ff = real(count(violation <= mu), dp)/real(size(violation), dp)
                mu = mu*(1.0_dp - ff/real(size(violation), dp))
            end if
        else
            ftrial = objective(trial)
            if (ieee_is_nan(ftrial)) error stop "jde_optimize: objective returned NaN"
            if (ftrial <= fpop(i)) then
                pop(:, i) = trial
                fpop(i) = ftrial
                hpop(:, i) = htrial
                violation(i) = trial_violation
                accepted = .true.
                ff = real(count(violation <= mu), dp)/real(size(violation), dp)
                mu = mu*(1.0_dp - ff/real(size(violation), dp))
            end if
        end if
    end subroutine constrained_selection_jde

    function continuation_rule(spread_value, tolerance, constrained, hpop, best) result(continue_search)
        real(dp), intent(in) :: spread_value, tolerance
        logical, intent(in) :: constrained
        real(dp), intent(in), optional :: hpop(:, :)
        integer, intent(in) :: best
        logical :: continue_search

        if (constrained) then
            if (.not. present(hpop)) error stop "continuation_rule: constraints missing"
            continue_search = spread_value >= tolerance .or. any(hpop(:, best) > 0.0_dp)
        else
            continue_search = spread_value >= tolerance
        end if
    end function continuation_rule

    function median_local(x) result(value)
        use deoptimr_utils, only: median_value
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = median_value(x)
    end function median_local
end module deoptimr_jde
