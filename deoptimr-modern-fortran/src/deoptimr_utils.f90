! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr_utils
    use deoptimr_kinds, only: dp
    use deoptimr_rng, only: random_uniform, random_integer, sample_without_replacement
    implicit none
    private

    public :: median_value, objective_spread, total_violation, transform_constraints
    public :: best_population_index, handle_bounds, make_population
    public :: initialize_adaptation, choose_distinct_indices, distance_squared
    public :: sort_indices_by_values, sort_indices_by_two_keys, nearest_distances
    public :: almost_equal, nearest_neighbor_index

contains

    function median_value(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        real(dp), allocatable :: work(:)
        integer :: i, j, n
        real(dp) :: key

        n = size(x)
        if (n == 0) error stop "median_value: empty input"
        allocate(work(n))
        work = x
        do i = 2, n
            key = work(i)
            j = i - 1
            do while (j >= 1)
                if (work(j) <= key) exit
                work(j + 1) = work(j)
                j = j - 1
            end do
            work(j + 1) = key
        end do
        if (mod(n, 2) == 1) then
            value = work((n + 1)/2)
        else
            value = 0.5_dp*(work(n/2) + work(n/2 + 1))
        end if
    end function median_value

    function objective_spread(cost, best, compare_to, scale) result(value)
        real(dp), intent(in) :: cost(:), best, scale
        character(len=*), intent(in) :: compare_to
        real(dp) :: value, reference

        if (scale <= 0.0_dp) error stop "objective_spread: scale must be positive"
        select case (trim(compare_to))
        case ("median")
            reference = median_value(cost)
        case ("max")
            reference = maxval(cost)
        case default
            error stop "objective_spread: compare_to must be median or max"
        end select
        value = (reference - best)/scale
    end function objective_spread

    pure function total_violation(values) result(value)
        real(dp), intent(in) :: values(:)
        real(dp) :: value
        value = sum(max(values, 0.0_dp))
    end function total_violation

    subroutine transform_constraints(raw, n_equalities, equality_tolerance, transformed)
        real(dp), intent(in) :: raw(:)
        integer, intent(in) :: n_equalities
        real(dp), intent(in), optional :: equality_tolerance(:)
        real(dp), intent(out) :: transformed(size(raw))
        integer :: i
        real(dp) :: eps

        if (n_equalities < 0 .or. n_equalities > size(raw)) then
            error stop "transform_constraints: invalid equality count"
        end if
        if (present(equality_tolerance)) then
            if (size(equality_tolerance) /= 1 .and. size(equality_tolerance) /= n_equalities) then
                error stop "transform_constraints: invalid tolerance length"
            end if
            if (any(equality_tolerance <= 0.0_dp)) then
                error stop "transform_constraints: tolerances must be positive"
            end if
        end if
        transformed = raw
        do i = 1, n_equalities
            eps = 1.0e-5_dp
            if (present(equality_tolerance)) then
                if (size(equality_tolerance) == 1) then
                    eps = equality_tolerance(1)
                else
                    eps = equality_tolerance(i)
                end if
            end if
            transformed(i) = abs(raw(i)) - eps
        end do
    end subroutine transform_constraints

    function best_population_index(cost, violation, mu, constrained) result(index_best)
        real(dp), intent(in) :: cost(:)
        real(dp), intent(in), optional :: violation(:), mu
        logical, intent(in) :: constrained
        integer :: index_best, i
        logical, allocatable :: feasible(:)
        real(dp) :: best_value

        if (.not. constrained) then
            index_best = minloc(cost, dim=1)
            return
        end if
        if (.not. present(violation) .or. .not. present(mu)) then
            error stop "best_population_index: missing constrained arguments"
        end if
        allocate(feasible(size(cost)))
        feasible = violation <= mu
        if (any(feasible)) then
            index_best = 0
            best_value = huge(1.0_dp)
            do i = 1, size(cost)
                if (feasible(i) .and. cost(i) < best_value) then
                    best_value = cost(i)
                    index_best = i
                end if
            end do
        else
            index_best = minloc(violation, dim=1)
        end if
    end function best_population_index

    subroutine handle_bounds(x, base, lower, upper)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: base(:), lower(:), upper(:)
        integer :: j

        if (size(x) /= size(base) .or. size(x) /= size(lower) .or. size(x) /= size(upper)) then
            error stop "handle_bounds: dimension mismatch"
        end if
        do j = 1, size(x)
            if (x(j) > upper(j)) x(j) = 0.5_dp*(upper(j) + base(j))
            if (x(j) < lower(j)) x(j) = 0.5_dp*(lower(j) + base(j))
        end do
    end subroutine handle_bounds

    subroutine make_population(lower, upper, n_random, population, initial_population)
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: n_random
        real(dp), allocatable, intent(out) :: population(:, :)
        real(dp), intent(in), optional :: initial_population(:, :)
        integer :: d, n_initial, i, j

        d = size(lower)
        if (size(upper) /= d) error stop "make_population: bound mismatch"
        if (n_random < 0) error stop "make_population: negative population size"
        n_initial = 0
        if (present(initial_population)) then
            if (size(initial_population, 1) /= d) error stop "make_population: initial row mismatch"
            if (any(initial_population < spread(lower, 2, size(initial_population, 2))) .or. &
                any(initial_population > spread(upper, 2, size(initial_population, 2)))) then
                error stop "make_population: initial point outside bounds"
            end if
            n_initial = size(initial_population, 2)
        end if
        allocate(population(d, n_random + n_initial))
        do j = 1, n_random
            do i = 1, d
                population(i, j) = random_uniform(lower(i), upper(i))
            end do
        end do
        if (n_initial > 0) population(:, n_random + 1:) = initial_population
    end subroutine make_population

    subroutine initialize_adaptation(d, np, f_lower, f_upper, cr_lower, cr_upper, &
            jitter_factor, use_jitter, f, cr, pf)
        integer, intent(in) :: d, np
        real(dp), intent(in) :: f_lower, f_upper, cr_lower, cr_upper, jitter_factor
        logical, intent(in) :: use_jitter
        real(dp), intent(out) :: f(:, :), cr(:), pf(:)
        integer :: i, j
        real(dp) :: scale

        if (size(f, 2) /= np .or. size(cr) /= np .or. size(pf) /= np) then
            error stop "initialize_adaptation: dimension mismatch"
        end if
        if (use_jitter) then
            if (size(f, 1) /= d) error stop "initialize_adaptation: jitter row mismatch"
            do j = 1, np
                scale = random_uniform(f_lower, f_upper)
                do i = 1, d
                    f(i, j) = scale*(1.0_dp + jitter_factor*random_uniform(-0.5_dp, 0.5_dp))
                end do
            end do
        else
            if (size(f, 1) /= 1) error stop "initialize_adaptation: dither row mismatch"
            do j = 1, np
                f(1, j) = random_uniform(f_lower, f_upper)
            end do
        end if
        do j = 1, np
            cr(j) = random_uniform(cr_lower, cr_upper)
            pf(j) = random_uniform()
        end do
    end subroutine initialize_adaptation

    subroutine choose_distinct_indices(np, excluded, selected)
        integer, intent(in) :: np, excluded
        integer, intent(out) :: selected(3)
        integer, allocatable :: candidates(:)
        integer :: i, k

        if (np < 4) error stop "choose_distinct_indices: population must be at least four"
        allocate(candidates(np - 1))
        k = 0
        do i = 1, np
            if (i /= excluded) then
                k = k + 1
                candidates(k) = i
            end if
        end do
        call sample_without_replacement(candidates, 3, selected)
    end subroutine choose_distinct_indices

    pure function distance_squared(x, y) result(value)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: value
        value = sum((x - y)**2)
    end function distance_squared

    function nearest_neighbor_index(x, population) result(index_nearest)
        real(dp), intent(in) :: x(:), population(:, :)
        integer :: index_nearest, j
        real(dp) :: best_distance, current

        index_nearest = 1
        best_distance = distance_squared(x, population(:, 1))
        do j = 2, size(population, 2)
            current = distance_squared(x, population(:, j))
            if (current < best_distance) then
                best_distance = current
                index_nearest = j
            end if
        end do
    end function nearest_neighbor_index

    subroutine sort_indices_by_values(values, indices)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: indices(size(values))
        integer :: i, j, key

        do i = 1, size(indices)
            indices(i) = i
        end do
        do i = 2, size(indices)
            key = indices(i)
            j = i - 1
            do while (j >= 1)
                if (values(indices(j)) <= values(key)) exit
                indices(j + 1) = indices(j)
                j = j - 1
            end do
            indices(j + 1) = key
        end do
    end subroutine sort_indices_by_values

    subroutine sort_indices_by_two_keys(primary, secondary, indices)
        real(dp), intent(in) :: primary(:), secondary(:)
        integer, intent(out) :: indices(size(primary))
        integer :: i, j, key

        if (size(secondary) /= size(primary)) error stop "sort_indices_by_two_keys: mismatch"
        do i = 1, size(indices)
            indices(i) = i
        end do
        do i = 2, size(indices)
            key = indices(i)
            j = i - 1
            do while (j >= 1)
                if (primary(indices(j)) < primary(key)) exit
                if (primary(indices(j)) > primary(key)) then
                    indices(j + 1) = indices(j)
                    j = j - 1
                    cycle
                end if
                if (secondary(indices(j)) <= secondary(key)) exit
                indices(j + 1) = indices(j)
                j = j - 1
            end do
            indices(j + 1) = key
        end do
    end subroutine sort_indices_by_two_keys

    subroutine nearest_distances(population, distances)
        real(dp), intent(in) :: population(:, :)
        real(dp), intent(out) :: distances(size(population, 2))
        integer :: i, j
        real(dp) :: current, best

        if (size(population, 2) < 2) error stop "nearest_distances: population too small"
        do i = 1, size(population, 2)
            best = huge(1.0_dp)
            do j = 1, size(population, 2)
                if (j == i) cycle
                current = distance_squared(population(:, i), population(:, j))
                if (current < best) best = current
            end do
            distances(i) = sqrt(best)
        end do
    end subroutine nearest_distances

    pure function almost_equal(a, b, tolerance) result(equal)
        real(dp), intent(in) :: a, b, tolerance
        logical :: equal
        real(dp) :: scale

        scale = max(1.0_dp, abs(a), abs(b))
        equal = abs(a - b) <= tolerance*scale
    end function almost_equal
end module deoptimr_utils
