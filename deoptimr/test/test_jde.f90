! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Portions derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This program is free software: you may redistribute it and/or modify it
! under the terms of GNU GPL version 2, or (at your option) any later version.
program test_jde
    use deoptimr, only: dp, seed_rng, jde_control, de_result, jde_optimize, spjde_optimize
    use deoptimr_utils, only: handle_bounds, transform_constraints, median_value
    use test_support, only: assert_true, assert_close, assert_all_close
    use test_benchmarks, only: sphere, quadratic, equality_constraint, square1
    use test_benchmarks, only: lower_one_constraint, aluffi
    implicit none

    call test_helpers()
    call test_jde_sphere()
    call test_jde_dither_max()
    call test_jde_equality()
    call test_jde_inequality()
    call test_spjde_aluffi()
    call test_spjde_equality()
    call test_spjde_inequality()
    call test_initial_population()
    write(*, '(a)') 'JDE and SPJDE tests passed.'

contains

    subroutine test_helpers()
        real(dp) :: x(2), base(2), lower(2), upper(2), raw(3), transformed(3), eps(2)
        x = [3.0_dp, -3.0_dp]
        base = [0.0_dp, 1.0_dp]
        lower = [-2.0_dp, -1.0_dp]
        upper = [2.0_dp, 3.0_dp]
        call handle_bounds(x, base, lower, upper)
        call assert_all_close(x, [1.0_dp, 0.0_dp], 1.0e-14_dp, 'bounce-back bounds')
        raw = [-0.2_dp, 0.03_dp, -1.0_dp]
        eps = [0.1_dp, 0.01_dp]
        call transform_constraints(raw, 2, eps, transformed)
        call assert_all_close(transformed, [0.1_dp, 0.02_dp, -1.0_dp], 1.0e-14_dp, &
            'equality constraint conversion')
        call assert_close(median_value([4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]), 2.5_dp, &
            1.0e-14_dp, 'median')
    end subroutine test_helpers

    subroutine test_jde_sphere()
        type(jde_control) :: control
        type(de_result) :: result

        call seed_rng(2345)
        control%population_size = 50
        control%max_iterations = 900
        control%tolerance = 1.0e-9_dp
        control%save_population = .true.
        call jde_optimize([-100.0_dp, -100.0_dp], [100.0_dp, 100.0_dp], sphere, result, control)
        call assert_true(result%value < 1.0e-10_dp, 'JDE sphere objective')
        call assert_true(maxval(abs(result%parameters)) < 2.0e-5_dp, 'JDE sphere minimizer')
        call assert_true(size(result%population, 2) == 50, 'JDE population details')
        call assert_true(result%convergence == 0, 'JDE convergence status')
    end subroutine test_jde_sphere


    subroutine test_jde_dither_max()
        type(jde_control) :: control
        type(de_result) :: result

        call seed_rng(1881)
        control%population_size = 45
        control%max_iterations = 1000
        control%tolerance = 1.0e-8_dp
        control%use_jitter = .false.
        control%compare_to = "max"
        call jde_optimize([-20.0_dp, -20.0_dp], [20.0_dp, 20.0_dp], sphere, result, control)
        call assert_true(result%value < 1.0e-8_dp, 'JDE dither-only objective')
        call assert_true(maxval(abs(result%parameters)) < 2.0e-4_dp, &
            'JDE maximum-spread convergence')
    end subroutine test_jde_dither_max

    subroutine test_jde_equality()
        type(jde_control) :: control
        type(de_result) :: result
        real(dp) :: eps(1)

        call seed_rng(7411)
        control%population_size = 70
        control%max_iterations = 1800
        control%tolerance = 1.0e-8_dp
        eps = 1.0e-5_dp
        call jde_optimize([-2.0_dp, -2.0_dp], [2.0_dp, 2.0_dp], quadratic, result, control, &
            equality_constraint, 1, 1, eps)
        call assert_true(abs(sum(result%parameters) - 1.0_dp) <= 2.0e-5_dp, &
            'JDE equality feasibility')
        call assert_all_close(result%parameters, [0.5_dp, 0.5_dp], 5.0e-3_dp, &
            'JDE equality minimizer')
        call assert_close(result%value, 0.5_dp, 5.0e-3_dp, 'JDE equality objective')
    end subroutine test_jde_equality

    subroutine test_jde_inequality()
        type(jde_control) :: control
        type(de_result) :: result

        call seed_rng(9001)
        control%population_size = 50
        control%max_iterations = 1200
        control%tolerance = 1.0e-9_dp
        call jde_optimize([-2.0_dp], [2.0_dp], square1, result, control, &
            lower_one_constraint, 1, 0)
        call assert_true(result%parameters(1) >= 1.0_dp - 2.0e-5_dp, &
            'JDE inequality feasibility')
        call assert_close(result%parameters(1), 1.0_dp, 2.0e-3_dp, 'JDE inequality minimizer')
        call assert_close(result%value, 1.0_dp, 4.0e-3_dp, 'JDE inequality objective')
    end subroutine test_jde_inequality

    subroutine test_spjde_aluffi()
        type(jde_control) :: control
        type(de_result) :: result

        call seed_rng(42)
        control%population_size = 40
        control%max_iterations = 1600
        control%tolerance = 1.0e-9_dp
        control%cr_lower = 0.0_dp
        control%cr_upper = 1.1_dp
        call spjde_optimize([-10.0_dp, -10.0_dp], [10.0_dp, 10.0_dp], aluffi, result, control)
        call assert_all_close(result%parameters, [-1.0466805_dp, 0.0_dp], 3.0e-3_dp, &
            'SPJDE Aluffi-Pentini minimizer')
        call assert_close(result%value, -0.352386_dp, 3.0e-3_dp, &
            'SPJDE Aluffi-Pentini objective')
    end subroutine test_spjde_aluffi


    subroutine test_spjde_equality()
        type(jde_control) :: control
        type(de_result) :: result
        real(dp) :: eps(1)

        call seed_rng(5521)
        control%population_size = 70
        control%max_iterations = 1800
        control%tolerance = 1.0e-8_dp
        control%cr_lower = 0.0_dp
        control%cr_upper = 1.1_dp
        eps = 1.0e-5_dp
        call spjde_optimize([-2.0_dp, -2.0_dp], [2.0_dp, 2.0_dp], quadratic, result, control, &
            equality_constraint, 1, 1, eps)
        call assert_true(abs(sum(result%parameters) - 1.0_dp) <= 2.0e-5_dp, &
            'SPJDE equality feasibility')
        call assert_all_close(result%parameters, [0.5_dp, 0.5_dp], 7.0e-3_dp, &
            'SPJDE equality minimizer')
    end subroutine test_spjde_equality


    subroutine test_spjde_inequality()
        type(jde_control) :: control
        type(de_result) :: result

        call seed_rng(7349)
        control%population_size = 60
        control%max_iterations = 1400
        control%tolerance = 1.0e-9_dp
        control%cr_upper = 1.1_dp
        call spjde_optimize([-2.0_dp], [2.0_dp], square1, result, control, &
            lower_one_constraint, 1, 0)
        call assert_true(result%parameters(1) >= 1.0_dp - 2.0e-5_dp, &
            'SPJDE inequality feasibility')
        call assert_close(result%parameters(1), 1.0_dp, 3.0e-3_dp, &
            'SPJDE inequality minimizer')
    end subroutine test_spjde_inequality

    subroutine test_initial_population()
        type(jde_control) :: control
        type(de_result) :: result
        real(dp) :: initial(2, 1)

        initial(:, 1) = [0.0_dp, 0.0_dp]
        call seed_rng(81)
        control%population_size = 3
        control%max_iterations = 0
        control%save_population = .true.
        call jde_optimize([-1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp], sphere, result, control, &
            initial_population=initial)
        call assert_true(size(result%population, 2) == 4, 'initial population appended')
        call assert_close(result%value, 0.0_dp, 1.0e-14_dp, 'initial optimum retained')
    end subroutine test_initial_population

end program test_jde
