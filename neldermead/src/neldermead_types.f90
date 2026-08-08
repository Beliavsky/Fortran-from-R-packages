! SPDX-License-Identifier: CECILL-2.0
! Derived from the R package neldermead 1.0-13 and its Scilab lineage.
! See LICENSE and UPSTREAM_PROVENANCE.md.

module neldermead_types
  use neldermead_kinds, only : dp
  implicit none
  private

  public :: nm_options, nm_result, nm_simplex
  public :: nm_objective, nm_constraints, nm_output_callback, nm_termination_callback
  public :: nm_method_variable, nm_method_fixed, nm_method_box
  public :: nm_simplex_given, nm_simplex_axes, nm_simplex_spendley
  public :: nm_simplex_pfeffer, nm_simplex_randbounds
  public :: nm_restart_oneill, nm_restart_kelley

  character(len=*), parameter :: nm_method_variable = 'variable'
  character(len=*), parameter :: nm_method_fixed    = 'fixed'
  character(len=*), parameter :: nm_method_box      = 'box'
  character(len=*), parameter :: nm_simplex_given   = 'given'
  character(len=*), parameter :: nm_simplex_axes    = 'axes'
  character(len=*), parameter :: nm_simplex_spendley= 'spendley'
  character(len=*), parameter :: nm_simplex_pfeffer = 'pfeffer'
  character(len=*), parameter :: nm_simplex_randbounds = 'randbounds'
  character(len=*), parameter :: nm_restart_oneill  = 'oneill'
  character(len=*), parameter :: nm_restart_kelley  = 'kelley'

  abstract interface
    function nm_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function nm_objective

    subroutine nm_constraints(x, c)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: c(:)
    end subroutine nm_constraints

    subroutine nm_output_callback(iteration, func_count, xbest, fbest, step, stop)
      import dp
      integer, intent(in) :: iteration, func_count
      real(dp), intent(in) :: xbest(:), fbest
      character(len=*), intent(in) :: step
      logical, intent(inout) :: stop
    end subroutine nm_output_callback

    subroutine nm_termination_callback(iteration, func_count, simplex, fv, terminate, status)
      import dp
      integer, intent(in) :: iteration, func_count
      real(dp), intent(in) :: simplex(:,:), fv(:)
      logical, intent(inout) :: terminate
      character(len=*), intent(inout) :: status
    end subroutine nm_termination_callback
  end interface

  type :: nm_options
    character(len=16) :: method = nm_method_variable
    character(len=16) :: simplex0_method = nm_simplex_axes
    character(len=16) :: restart_simplex_method = 'oriented'
    character(len=16) :: restart_detection = nm_restart_oneill
    integer :: max_iter = 1000
    integer :: max_fun_evals = 2000
    real(dp) :: rho = 1.0_dp
    real(dp) :: chi = 2.0_dp
    real(dp) :: gamma = 0.5_dp
    real(dp) :: sigma = 0.5_dp
    logical :: greedy = .false.
    real(dp) :: simplex0_length = 1.0_dp
    real(dp) :: simplex0_delta_usual = 0.05_dp
    real(dp) :: simplex0_delta_zero = 0.0075_dp
    integer :: box_npoints = 0
    real(dp) :: box_reflect = 1.3_dp
    real(dp) :: box_ineq_scaling = 0.5_dp
    real(dp) :: box_bounds_alpha = 1.0e-6_dp
    real(dp) :: gui_alpha_min = 1.0e-5_dp
    logical :: tol_x_method = .true.
    logical :: tol_fun_method = .true.
    real(dp) :: tol_x_absolute = 0.0_dp
    real(dp) :: tol_x_relative = epsilon(1.0_dp)
    real(dp) :: tol_fun_absolute = 0.0_dp
    real(dp) :: tol_fun_relative = epsilon(1.0_dp)
    logical :: tol_simplex_size_method = .true.
    real(dp) :: tol_simplex_size_absolute = 0.0_dp
    real(dp) :: tol_simplex_size_relative = epsilon(1.0_dp)
    logical :: tol_size_delta_f_method = .false.
    real(dp) :: tol_delta_f = epsilon(1.0_dp)
    logical :: tol_f_std_method = .false.
    real(dp) :: tol_f_std = 0.0_dp
    logical :: tol_variance_flag = .false.
    real(dp) :: tol_absolute_variance = 0.0_dp
    real(dp) :: tol_relative_variance = epsilon(1.0_dp)
    logical :: box_termination = .false.
    real(dp) :: box_tol_f = 1.0e-5_dp
    integer :: box_nb_match = 5
    logical :: kelley_stagnation_flag = .false.
    logical :: kelley_normalization_flag = .true.
    real(dp) :: kelley_stagnation_alpha0 = 1.0e-4_dp
    logical :: restart_flag = .false.
    integer :: restart_max = 3
    real(dp) :: restart_eps = epsilon(1.0_dp)
    real(dp) :: restart_step = 1.0_dp
    logical :: store_history = .false.
    logical :: verbose = .false.
    integer :: seed = 12345
  end type nm_options

  type :: nm_simplex
    real(dp), allocatable :: x(:,:)   ! (n, nvertices)
    real(dp), allocatable :: f(:)
  end type nm_simplex

  type :: nm_result
    real(dp), allocatable :: x(:)
    real(dp) :: f = huge(1.0_dp)
    integer :: iterations = 0
    integer :: func_count = 0
    integer :: restart_count = 0
    logical :: converged = .false.
    character(len=32) :: status = 'not_started'
    character(len=48) :: algorithm = ''
    type(nm_simplex) :: simplex
    real(dp), allocatable :: history_x(:,:)
    real(dp), allocatable :: history_f(:)
    integer :: history_count = 0
  end type nm_result

end module neldermead_types
