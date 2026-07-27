! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr_types
    use deoptimr_kinds, only: dp
    implicit none
    private

    public :: jde_control, ncde_control, de_result, ncde_result

    type :: jde_control
        integer :: population_size = 0
        real(dp) :: f_lower = 0.1_dp
        real(dp) :: f_upper = 1.0_dp
        real(dp) :: cr_lower = 0.0_dp
        real(dp) :: cr_upper = 1.0_dp
        real(dp) :: tau_f = 0.1_dp
        real(dp) :: tau_cr = 0.1_dp
        real(dp) :: tau_pf = 0.1_dp
        real(dp) :: jitter_factor = 0.001_dp
        logical :: use_jitter = .true.
        real(dp) :: tolerance = 1.0e-15_dp
        integer :: max_iterations = 0
        real(dp) :: objective_scale = 1.0_dp
        character(len=6) :: compare_to = "median"
        logical :: trace = .false.
        integer :: trace_interval = 1
        logical :: save_population = .true.
    end type jde_control

    type :: ncde_control
        integer :: population_size = 100
        real(dp) :: critical_tolerance = 1.0e-5_dp
        real(dp) :: niche_radius = -1.0_dp
        integer :: archive_size = 100
        logical :: reinitialize_archive_neighbors = .true.
        real(dp) :: f_lower = 0.1_dp
        real(dp) :: f_upper = 1.0_dp
        real(dp) :: cr_lower = 0.0_dp
        real(dp) :: cr_upper = 1.1_dp
        integer :: neighbor_lower = 0
        integer :: neighbor_upper = 0
        real(dp) :: tau_f = 0.1_dp
        real(dp) :: tau_cr = 0.1_dp
        real(dp) :: tau_pf = 0.1_dp
        real(dp) :: tau_neighbors = 0.1_dp
        real(dp) :: jitter_factor = 0.001_dp
        logical :: use_jitter = .true.
        integer :: max_iterations = 2000
        logical :: trace = .false.
        integer :: trace_interval = 1
    end type ncde_control

    type :: de_result
        real(dp), allocatable :: parameters(:)
        real(dp) :: value = huge(1.0_dp)
        integer :: iterations = 0
        integer :: convergence = 1
        real(dp), allocatable :: population(:, :)
        real(dp), allocatable :: population_cost(:)
        real(dp), allocatable :: population_constraints(:, :)
        real(dp), allocatable :: total_violation(:)
    end type de_result

    type :: ncde_result
        integer :: iterations = 0
        real(dp), allocatable :: solution_archive(:, :)
        real(dp), allocatable :: objective_archive(:)
        real(dp), allocatable :: constraint_archive(:, :)
        real(dp), allocatable :: solution_population(:, :)
        real(dp), allocatable :: objective_population(:)
        real(dp), allocatable :: constraint_population(:, :)
        real(dp), allocatable :: total_violation_population(:)
        real(dp) :: final_niche_radius = huge(1.0_dp)
    end type ncde_result
end module deoptimr_types
