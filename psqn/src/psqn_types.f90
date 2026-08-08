! SPDX-License-Identifier: Apache-2.0
module psqn_types
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  private

  integer, parameter, public :: dp = real64

  integer, parameter, public :: psqn_max_it_reached = -1
  integer, parameter, public :: psqn_conjugate_gradient_failed = -2
  integer, parameter, public :: psqn_line_search_failed = -3
  integer, parameter, public :: psqn_user_interrupt = -4
  integer, parameter, public :: psqn_converged = 0

  integer, parameter, public :: psqn_pre_none = 0
  integer, parameter, public :: psqn_pre_diag = 1
  integer, parameter, public :: psqn_pre_cholesky = 2
  integer, parameter, public :: psqn_pre_block = 3

  type, public :: psqn_options
    real(dp) :: rel_eps = 1.0e-8_dp
    integer :: max_it = 100
    real(dp) :: c1 = 1.0e-4_dp
    real(dp) :: c2 = 0.9_dp
    logical :: use_bfgs = .true.
    integer :: trace = 0
    real(dp) :: cg_tol = 0.5_dp
    logical :: strong_wolfe = .true.
    integer :: max_cg = 0
    integer :: pre_method = psqn_pre_diag
    real(dp) :: gr_tol = -1.0_dp
  end type psqn_options

  type, public :: psqn_bfgs_options
    real(dp) :: rel_eps = 1.0e-8_dp
    integer :: max_it = 100
    real(dp) :: c1 = 1.0e-4_dp
    real(dp) :: c2 = 0.9_dp
    integer :: trace = 0
    real(dp) :: gr_tol = -1.0_dp
    real(dp) :: abs_eps = -1.0_dp
  end type psqn_bfgs_options

  type, public :: psqn_auglag_options
    real(dp) :: penalty_start = 1.0_dp
    integer :: max_it_outer = 100
    real(dp) :: violations_norm_thresh = 1.0e-6_dp
    real(dp) :: tau = 1.5_dp
  end type psqn_auglag_options

  type, public :: psqn_info
    real(dp) :: value = huge(1.0_dp)
    integer :: info = psqn_max_it_reached
    integer :: n_eval = 0
    integer :: n_grad = 0
    integer :: n_cg = 0
  end type psqn_info

  type, public :: psqn_auglag_info
    real(dp) :: value = huge(1.0_dp)
    integer :: info = psqn_max_it_reached
    integer :: n_eval = 0
    integer :: n_grad = 0
    integer :: n_cg = 0
    integer :: n_aug_lagrang = 0
    real(dp) :: penalty = 0.0_dp
  end type psqn_auglag_info

  type, public :: psqn_element_spec
    integer, allocatable :: idx(:)
  end type psqn_element_spec

  abstract interface
    subroutine psqn_element_eval(i, x, f, g, comp_grad)
      import dp
      integer, intent(in) :: i
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      logical, intent(in) :: comp_grad
    end subroutine psqn_element_eval

    subroutine psqn_objective_eval(x, f, g, comp_grad)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      logical, intent(in) :: comp_grad
    end subroutine psqn_objective_eval
  end interface

  public :: psqn_element_eval, psqn_objective_eval

end module psqn_types
