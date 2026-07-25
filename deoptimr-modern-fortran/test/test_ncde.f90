! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Portions derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This program is free software: you may redistribute it and/or modify it
! under the terms of GNU GPL version 2, or (at your option) any later version.
program test_ncde
    use deoptimr, only: dp, seed_rng, ncde_control, ncde_result, ncde_optimize
    use test_support, only: assert_true
    use test_benchmarks, only: becker_lago, square1, outside_unit_constraint
    use test_benchmarks, only: unit_equality_constraint
    implicit none

    call test_becker_lago()
    call test_constrained_two_minima()
    call test_equality_two_minima()
    call test_automatic_radius()
    write(*, '(a)') 'NCDE tests passed.'

contains

    subroutine test_becker_lago()
        type(ncde_control) :: control
        type(ncde_result) :: result
        logical :: found(4)
        integer :: j, quadrant

        call seed_rng(2345)
        control%population_size = 120
        control%niche_radius = 4.0_dp
        control%archive_size = 20
        control%max_iterations = 350
        control%neighbor_lower = 6
        control%neighbor_upper = 24
        call ncde_optimize([-10.0_dp, -10.0_dp], [10.0_dp, 10.0_dp], becker_lago, result, control)
        found = .false.
        do j = 1, size(result%objective_archive)
            if (result%objective_archive(j) > 2.0e-4_dp) cycle
            quadrant = merge(1, 0, result%solution_archive(1, j) > 0.0_dp) + &
                2*merge(1, 0, result%solution_archive(2, j) > 0.0_dp) + 1
            found(quadrant) = maxval(abs(abs(result%solution_archive(:, j)) - 5.0_dp)) < 2.0e-2_dp
        end do
        call assert_true(all(found), 'NCDE found all four Becker-Lago minima')
        call assert_true(size(result%solution_population, 2) == 120, 'NCDE final population')
    end subroutine test_becker_lago

    subroutine test_constrained_two_minima()
        type(ncde_control) :: control
        type(ncde_result) :: result
        logical :: found_negative, found_positive
        integer :: j

        call seed_rng(9182)
        control%population_size = 100
        control%niche_radius = 0.6_dp
        control%archive_size = 10
        control%max_iterations = 450
        control%neighbor_lower = 5
        control%neighbor_upper = 20
        call ncde_optimize([-2.0_dp], [2.0_dp], square1, result, control, &
            outside_unit_constraint, 1, 0)
        found_negative = .false.
        found_positive = .false.
        do j = 1, size(result%objective_archive)
            if (abs(abs(result%solution_archive(1, j)) - 1.0_dp) < 2.0e-2_dp .and. &
                result%objective_archive(j) < 1.05_dp) then
                if (result%solution_archive(1, j) < 0.0_dp) found_negative = .true.
                if (result%solution_archive(1, j) > 0.0_dp) found_positive = .true.
            end if
        end do
        call assert_true(found_negative .and. found_positive, &
            'NCDE found both constrained minima')
        call assert_true(all(result%constraint_archive <= 0.0_dp), &
            'NCDE archive is strictly feasible')
    end subroutine test_constrained_two_minima


    subroutine test_equality_two_minima()
        type(ncde_control) :: control
        type(ncde_result) :: result
        real(dp) :: eps(1)
        logical :: found_negative, found_positive
        integer :: j

        call seed_rng(12077)
        control%population_size = 120
        control%niche_radius = 0.6_dp
        control%archive_size = 10
        control%max_iterations = 650
        control%neighbor_lower = 6
        control%neighbor_upper = 24
        eps = 1.0e-4_dp
        call ncde_optimize([-2.0_dp], [2.0_dp], square1, result, control, &
            unit_equality_constraint, 1, 1, eps)
        found_negative = .false.
        found_positive = .false.
        do j = 1, size(result%objective_archive)
            if (abs(abs(result%solution_archive(1, j)) - 1.0_dp) < 2.0e-3_dp) then
                if (result%solution_archive(1, j) < 0.0_dp) found_negative = .true.
                if (result%solution_archive(1, j) > 0.0_dp) found_positive = .true.
            end if
        end do
        call assert_true(found_negative .and. found_positive, &
            'NCDE found both equality-constrained minima')
        call assert_true(all(result%constraint_archive <= 0.0_dp), &
            'NCDE equality archive is feasible')
    end subroutine test_equality_two_minima

    subroutine test_automatic_radius()
        type(ncde_control) :: control
        type(ncde_result) :: result

        call seed_rng(301)
        control%population_size = 50
        control%niche_radius = -1.0_dp
        control%archive_size = 5
        control%max_iterations = 10
        control%reinitialize_archive_neighbors = .false.
        control%neighbor_lower = 3
        control%neighbor_upper = 10
        call ncde_optimize([-4.0_dp, -4.0_dp], [4.0_dp, 4.0_dp], becker_lago, result, control)
        call assert_true(result%final_niche_radius > 0.0_dp, 'automatic niche radius positive')
        call assert_true(result%final_niche_radius < huge(1.0_dp), 'automatic niche radius finite')
    end subroutine test_automatic_radius

end program test_ncde
