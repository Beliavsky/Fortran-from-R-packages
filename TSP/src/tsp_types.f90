! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp_types
    use tsp_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: tsp_identity = 1
    integer, parameter, public :: tsp_random = 2
    integer, parameter, public :: tsp_nearest_insertion = 3
    integer, parameter, public :: tsp_farthest_insertion = 4
    integer, parameter, public :: tsp_cheapest_insertion = 5
    integer, parameter, public :: tsp_arbitrary_insertion = 6
    integer, parameter, public :: tsp_nn = 7
    integer, parameter, public :: tsp_repetitive_nn = 8
    integer, parameter, public :: tsp_two_opt_method = 9
    integer, parameter, public :: tsp_sa_method = 10

    integer, parameter, public :: sa_reversal = 1
    integer, parameter, public :: sa_swap = 2
    integer, parameter, public :: sa_mixed = 3

    type, public :: tsp_control
        integer :: start = 0
        integer :: rep = 1
        logical :: two_opt = .false.
        integer :: two_opt_repetitions = 1
        integer :: sa_move = sa_reversal
        real(dp) :: temp = -1.0_dp
        integer :: tmax = 10
        integer :: maxit = 10000
        integer, allocatable :: tour(:)
    end type tsp_control

    type, public :: tsp_tour
        integer, allocatable :: order(:)
        real(dp) :: length = 0.0_dp
        character(len=48) :: method = "unknown"
    end type tsp_tour

    type, public :: tsp_path
        integer, allocatable :: city(:)
    end type tsp_path

    type, public :: tsp_path_collection
        type(tsp_path), allocatable :: path(:)
    end type tsp_path_collection

end module tsp_types
