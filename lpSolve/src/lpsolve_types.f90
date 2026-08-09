! SPDX-License-Identifier: LGPL-2.0-only
module lpsolve_types
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    integer, parameter, public :: LP_MIN = 0
    integer, parameter, public :: LP_MAX = 1

    integer, parameter, public :: LP_LE = 1
    integer, parameter, public :: LP_GE = 2
    integer, parameter, public :: LP_EQ = 3

    integer, parameter, public :: LP_OPTIMAL = 0
    integer, parameter, public :: LP_SUBOPTIMAL = 1
    integer, parameter, public :: LP_INFEASIBLE = 2
    integer, parameter, public :: LP_UNBOUNDED = 3
    integer, parameter, public :: LP_NUMFAILURE = 5
    integer, parameter, public :: LP_TIMEOUT = 7

    real(dp), parameter, public :: LP_INFINITY = huge(1.0_dp) / 100.0_dp

    type, public :: lp_control
        real(dp) :: feasibility_tol = 1.0e-9_dp
        real(dp) :: optimality_tol = 1.0e-10_dp
        real(dp) :: integrality_tol = 1.0e-8_dp
        integer :: max_simplex_iter = 200000
        integer :: max_nodes = 100000
        real(dp) :: timeout_seconds = 0.0_dp
        logical :: scale_rows = .true.
        logical :: bland_rule = .true.
        integer :: num_binary_solutions = 1
    end type lp_control

    type, public :: lp_result
        integer :: status = LP_NUMFAILURE
        real(dp) :: objective = 0.0_dp
        real(dp), allocatable :: solution(:)
        real(dp), allocatable :: duals(:)
        real(dp), allocatable :: reduced_costs(:)
        real(dp), allocatable :: solutions(:,:)
        real(dp), allocatable :: solution_objectives(:)
        integer :: solution_count = 0
        integer :: simplex_iterations = 0
        integer :: nodes = 0
        logical :: sensitivity_ranges_available = .false.
    end type lp_result

    type, public :: sparse_constraints
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: row(:)
        integer, allocatable :: col(:)
        real(dp), allocatable :: val(:)
    end type sparse_constraints

end module lpsolve_types
