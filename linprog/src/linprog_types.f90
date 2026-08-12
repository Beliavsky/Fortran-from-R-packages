! SPDX-License-Identifier: GPL-2.0-or-later
module linprog_types
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    integer, parameter, public :: LINPROG_LE = -1
    integer, parameter, public :: LINPROG_EQ = 0
    integer, parameter, public :: LINPROG_GE = 1

    integer, parameter, public :: LINPROG_SUCCESS = 0
    integer, parameter, public :: LINPROG_LPSOLVE_FAILED = 1
    integer, parameter, public :: LINPROG_DUAL_FAILED = 2
    integer, parameter, public :: LINPROG_CONSTRAINT_VIOLATION = 3
    integer, parameter, public :: LINPROG_PHASE1_MAXITER = 4
    integer, parameter, public :: LINPROG_PHASE2_MAXITER = 5

    type, public :: linprog_control
        logical :: maximum = .false.
        integer :: maxiter = 1000
        real(dp) :: zero = 1.0e-9_dp
        real(dp) :: tol = 1.0e-6_dp
        real(dp) :: dualtol = 1.0e-6_dp
        logical :: use_lpsolve = .false.
        logical :: solve_dual = .false.
        integer :: verbose = 0
    end type linprog_control

    type, public :: linprog_result
        integer :: status = LINPROG_SUCCESS
        integer :: lp_status = 0
        integer :: dual_status = 0
        real(dp) :: opt = 0.0_dp
        integer :: iter1 = 0
        integer :: iter2 = 0
        logical :: maximum = .false.
        logical :: used_lpsolve = .false.
        logical :: solved_dual = .false.
        integer :: maxiter = 1000
        real(dp), allocatable :: solution(:)
        integer, allocatable :: basis(:)
        real(dp), allocatable :: basic_values(:)
        real(dp), allocatable :: allvar_opt(:)
        real(dp), allocatable :: allvar_cvec(:)
        real(dp), allocatable :: allvar_min_c(:)
        real(dp), allocatable :: allvar_max_c(:)
        real(dp), allocatable :: allvar_marg(:)
        real(dp), allocatable :: allvar_marg_reg(:)
        real(dp), allocatable :: con_actual(:)
        real(dp), allocatable :: con_bvec(:)
        real(dp), allocatable :: con_free(:)
        real(dp), allocatable :: con_dual(:)
        real(dp), allocatable :: con_dual_reg(:)
        real(dp), allocatable :: con_dual_primal(:)
        integer, allocatable :: con_dir(:)
        real(dp), allocatable :: tableau(:,:)
    end type linprog_result

    type, public :: mps_model
        character(len=:), allocatable :: name
        real(dp), allocatable :: cvec(:)
        real(dp), allocatable :: bvec(:)
        real(dp), allocatable :: amat(:,:)
        character(len=:), allocatable :: var_names(:)
        character(len=:), allocatable :: con_names(:)
        logical :: has_result = .false.
        type(linprog_result) :: result
    end type mps_model

    public :: linprog_nan

contains

    pure real(dp) function linprog_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function linprog_nan

end module linprog_types
