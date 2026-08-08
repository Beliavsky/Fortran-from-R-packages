! SPDX-License-Identifier: MIT
module optimflex_types
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  integer, parameter, public :: diff_forward = 1
  integer, parameter, public :: diff_central = 2
  integer, parameter, public :: diff_richardson = 3

  type, public :: optim_control
    logical :: use_abs_f = .false.
    logical :: use_rel_f = .false.
    logical :: use_abs_x = .false.
    logical :: use_rel_x = .true.
    logical :: use_grad = .true.
    logical :: use_posdef = .true.
    logical :: use_pred_f = .false.
    logical :: use_pred_f_avg = .false.
    logical :: use_damped = .true.
    integer :: max_iter = 1000
    integer :: diff_method = diff_forward
    integer :: ls_max_steps = 30
    integer :: zoom_max_steps = 25
    integer :: memory = 5
    integer :: ridge_max_tries = 8
    integer :: gonfle_max = 10
    real(dp) :: tol_abs_f = 1.0e-6_dp
    real(dp) :: tol_rel_f = 1.0e-6_dp
    real(dp) :: tol_abs_x = 1.0e-6_dp
    real(dp) :: tol_rel_x = 1.0e-6_dp
    real(dp) :: tol_grad = 1.0e-4_dp
    real(dp) :: tol_pred_f = 1.0e-4_dp
    real(dp) :: tol_pred_f_avg = 1.0e-4_dp
    real(dp) :: wolfe_c1 = 1.0e-4_dp
    real(dp) :: wolfe_c2 = 0.9_dp
    real(dp) :: ls_alpha0 = 1.0_dp
    real(dp) :: ls_shrink = 0.5_dp
    real(dp) :: curvature_eps = 1.0e-12_dp
    real(dp) :: hinv_init_diag = 1.0_dp
    real(dp) :: h_init_diag = 1.0_dp
    real(dp) :: damp_phi = 0.2_dp
    real(dp) :: bound_eps = 1.0e-10_dp
    real(dp) :: ridge_offset = 1.0e-4_dp
    real(dp) :: ridge_mult = 10.0_dp
    real(dp) :: initial_delta = 1.0_dp
    real(dp) :: delta_max = 100.0_dp
    real(dp) :: rho_accept = 0.1_dp
    real(dp) :: rho_expand = 0.75_dp
    real(dp) :: delta_shrink = 0.25_dp
    real(dp) :: delta_expand = 2.0_dp
    real(dp) :: dd_bias = 0.8_dp
    real(dp) :: ridge_eps = 1.0e-9_dp
    real(dp) :: da_init = 1.0e-2_dp
    real(dp) :: ga_init = 1.0e-2_dp
    real(dp) :: da_factor = 5.0_dp
    real(dp) :: da_min = 1.0e-7_dp
    real(dp) :: ls_min_step = 1.0e-14_dp
    character(len=16) :: hessian_update = 'bfgs'
  end type optim_control

  type, public :: optim_result
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: approx_hessian(:,:)
    real(dp), allocatable :: approx_hinv(:,:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: cpu_time = 0.0_dp
    real(dp) :: elapsed_time = 0.0_dp
    real(dp) :: max_grad = huge(1.0_dp)
    real(dp) :: pred_dec = huge(1.0_dp)
    real(dp) :: pred_dec_avg = huge(1.0_dp)
    logical :: converged = .false.
    logical :: hess_is_pd = .false.
    integer :: iter = 0
    character(len=64) :: status = 'not_started'
  end type optim_result

  abstract interface
    function objective_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_fn

    subroutine gradient_fn(x, g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine gradient_fn

    subroutine hessian_fn(x, h)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: h(:,:)
    end subroutine hessian_fn

    function residual_fn(x) result(r)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: r(:)
    end function residual_fn

    function jacobian_fn(x) result(j)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:,:)
    end function jacobian_fn
  end interface

  public :: objective_fn, gradient_fn, hessian_fn, residual_fn, jacobian_fn
  public :: bfgs_default_control, lbfgsb_default_control
  public :: newton_default_control, modified_newton_default_control
  public :: gauss_newton_default_control, lm_default_control
  public :: dogleg_default_control, double_dogleg_default_control

contains

  function bfgs_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 10000
    c%curvature_eps = 1.0e-12_dp
    c%hinv_init_diag = 1.0_dp
  end function bfgs_default_control

  function lbfgsb_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 10000
    c%memory = 5
    c%curvature_eps = 1.0e-10_dp
    c%bound_eps = 1.0e-10_dp
  end function lbfgsb_default_control

  function newton_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 1000
  end function newton_default_control

  function modified_newton_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 1000
    c%ridge_offset = 1.0e-4_dp
  end function modified_newton_default_control

  function gauss_newton_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 1000
    c%ridge_offset = 1.0e-4_dp
    c%ridge_mult = 10.0_dp
    c%ridge_max_tries = 8
  end function gauss_newton_default_control

  function lm_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 10000
    c%da_init = 1.0e-2_dp
    c%ga_init = 1.0e-2_dp
    c%da_factor = 5.0_dp
    c%da_min = 1.0e-7_dp
    c%gonfle_max = 10
    c%ls_max_steps = 20
    c%ls_min_step = 1.0e-14_dp
  end function lm_default_control

  function dogleg_default_control() result(c)
    type(optim_control) :: c
    c%max_iter = 10000
    c%initial_delta = 1.0_dp
    c%delta_max = 100.0_dp
  end function dogleg_default_control

  function double_dogleg_default_control() result(c)
    type(optim_control) :: c
    c = dogleg_default_control()
    c%dd_bias = 0.8_dp
    c%ridge_eps = 1.0e-9_dp
  end function double_dogleg_default_control

end module optimflex_types
