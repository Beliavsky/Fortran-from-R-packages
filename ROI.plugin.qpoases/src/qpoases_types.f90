! SPDX-License-Identifier: GPL-3.0-only
module qpoases_types
    use qpoases_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: hst_zero = 0
    integer, parameter, public :: hst_identity = 1
    integer, parameter, public :: hst_posdef = 2
    integer, parameter, public :: hst_posdef_nullspace = 3
    integer, parameter, public :: hst_semidef = 4
    integer, parameter, public :: hst_indef = 5
    integer, parameter, public :: hst_unknown = 6

    integer, parameter, public :: successful_return = 0
    integer, parameter, public :: ret_invalid_arguments = 3
    integer, parameter, public :: ret_init_failed = 33
    integer, parameter, public :: ret_init_failed_infeasibility = 37
    integer, parameter, public :: ret_init_failed_unboundedness = 38
    integer, parameter, public :: ret_qp_unbounded = 45
    integer, parameter, public :: ret_qp_infeasible = 46
    integer, parameter, public :: ret_qp_not_solved = 47
    integer, parameter, public :: ret_qp_solved = 48
    integer, parameter, public :: ret_hotstart_failed = 51
    integer, parameter, public :: ret_max_nwsr_reached = 63
    integer, parameter, public :: ret_hessian_not_spd = 98

    type, public :: qpoases_options
        integer :: print_level = 1
        logical :: enable_ramping = .true.
        logical :: enable_far_bounds = .true.
        logical :: enable_flipping_bounds = .true.
        logical :: enable_regularisation = .true.
        logical :: enable_full_li_tests = .false.
        logical :: enable_nzc_tests = .true.
        integer :: enable_drift_correction = 1
        integer :: enable_cholesky_refactorisation = 0
        logical :: enable_equalities = .false.
        real(dp) :: termination_tolerance = 5.0e6_dp * epsilon(1.0_dp)
        real(dp) :: bound_tolerance = 1.0e6_dp * epsilon(1.0_dp)
        real(dp) :: bound_relaxation = 1.0e4_dp
        real(dp) :: eps_num = -1.0e3_dp * epsilon(1.0_dp)
        real(dp) :: eps_den = 1.0e3_dp * epsilon(1.0_dp)
        real(dp) :: max_primal_jump = 1.0e8_dp
        real(dp) :: max_dual_jump = 1.0e8_dp
        real(dp) :: initial_ramping = 0.5_dp
        real(dp) :: final_ramping = 1.0_dp
        real(dp) :: initial_far_bounds = 1.0e6_dp
        real(dp) :: grow_far_bounds = 1.0e3_dp
        real(dp) :: rcond_s_min = 1.0e-14_dp
        real(dp) :: eps_flipping = 1.0e3_dp * epsilon(1.0_dp)
        real(dp) :: eps_regularisation = 1.0e3_dp * epsilon(1.0_dp)
        real(dp) :: eps_iter_ref = 1.0e2_dp * epsilon(1.0_dp)
        real(dp) :: eps_li_tests = 1.0e5_dp * epsilon(1.0_dp)
        real(dp) :: eps_nzc_tests = 3.0e3_dp * epsilon(1.0_dp)
        integer :: num_regularisation_steps = 1
        integer :: num_refinement_steps = 1
        integer :: drop_bound_priority = 1
        integer :: drop_eq_con_priority = 1
        integer :: drop_ineq_con_priority = 1
        logical :: enable_inertia_correction = .true.
        logical :: enable_drop_infeasibles = .false.
        integer :: initial_status_bounds = -1
    end type qpoases_options

    type, public :: qpoases_result
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp) :: objval = huge(1.0_dp)
        integer :: status = ret_qp_not_solved
        integer :: nwsr = 0
        logical :: initialised = .false.
        logical :: solved = .false.
        logical :: infeasible = .false.
        logical :: unbounded = .false.
        integer :: n_free = 0
        integer :: n_fixed = 0
        integer :: n_active_constraints = 0
        integer :: n_inactive_constraints = 0
    end type qpoases_result
end module qpoases_types
