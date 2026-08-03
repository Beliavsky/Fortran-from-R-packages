! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_types
  use maxlik_kinds, only: dp
  use maxlik_status, only: MAXLIK_INVALID_INPUT
  implicit none
  private

  abstract interface
    subroutine maxlik_objective(x, value, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      integer, intent(out) :: status
    end subroutine maxlik_objective

    subroutine maxlik_gradient(x, gradient, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      integer, intent(out) :: status
    end subroutine maxlik_gradient

    subroutine maxlik_hessian(x, hessian, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hessian(:, :)
      integer, intent(out) :: status
    end subroutine maxlik_hessian

    subroutine maxlik_scores(x, scores, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: scores(:, :)
      integer, intent(out) :: status
    end subroutine maxlik_scores
  end interface

  type, public :: maxlik_problem
    integer :: npar = 0
    integer :: nobs = 0
    procedure(maxlik_objective), pointer, nopass :: objective => null()
    procedure(maxlik_gradient), pointer, nopass :: gradient => null()
    procedure(maxlik_hessian), pointer, nopass :: hessian => null()
    procedure(maxlik_scores), pointer, nopass :: scores => null()
    logical, allocatable :: active(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    real(dp), allocatable :: eq_a(:, :)
    real(dp), allocatable :: eq_b(:)
    real(dp), allocatable :: ineq_a(:, :)
    real(dp), allocatable :: ineq_b(:)
    real(dp) :: penalty_rho = 0.0_dp
  end type maxlik_problem

  type, public :: maxlik_control
    real(dp) :: tol = 1.0e-8_dp
    real(dp) :: reltol = sqrt(epsilon(1.0_dp))
    real(dp) :: gradtol = 1.0e-6_dp
    real(dp) :: steptol = 1.0e-10_dp
    real(dp) :: lambdatol = 1.0e-6_dp
    real(dp) :: qrtol = 1.0e-10_dp
    character(len=16) :: qac = 'stephalving'
    real(dp) :: marquardt_lambda0 = 1.0e-2_dp
    real(dp) :: marquardt_lambda_step = 2.0_dp
    real(dp) :: marquardt_max_lambda = 1.0e12_dp
    real(dp) :: nm_alpha = 1.0_dp
    real(dp) :: nm_beta = 0.5_dp
    real(dp) :: nm_gamma = 2.0_dp
    real(dp) :: sann_temp = 10.0_dp
    integer :: sann_tmax = 10
    integer :: random_seed = 123
    real(dp) :: sga_momentum = 0.0_dp
    real(dp) :: adam_momentum1 = 0.9_dp
    real(dp) :: adam_momentum2 = 0.999_dp
    real(dp) :: learning_rate = 0.1_dp
    integer :: batch_size = 0
    real(dp) :: gradient_clip = 0.0_dp
    integer :: patience = 0
    integer :: patience_step = 1
    integer :: iterlim = 150
    integer :: print_level = 0
    logical :: store_values = .false.
    logical :: store_parameters = .false.
    logical :: final_hessian = .true.
    logical :: use_central_differences = .true.
    integer :: constraint_max_outer = 20
    real(dp) :: constraint_tol = sqrt(epsilon(1.0_dp))
    real(dp) :: constraint_rho0 = 1.0_dp
    real(dp) :: constraint_rho_factor = 10.0_dp
  end type maxlik_control

  type, public :: maxlik_result
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: std_error(:)
    real(dp), allocatable :: gradient_obs(:, :)
    real(dp), allocatable :: stored_values(:)
    real(dp), allocatable :: stored_parameters(:, :)
    logical, allocatable :: active(:)
    real(dp) :: maximum = -huge(1.0_dp)
    real(dp) :: condition_number = huge(1.0_dp)
    real(dp) :: constraint_violation = 0.0_dp
    integer :: code = MAXLIK_INVALID_INPUT
    integer :: iterations = 0
    integer :: outer_iterations = 0
    integer :: function_count = 0
    integer :: gradient_count = 0
    integer :: hessian_count = 0
    logical :: converged = .false.
    character(len=40) :: method = ''
    character(len=96) :: message = 'not run'
  end type maxlik_result

  type, public :: derivative_comparison
    real(dp), allocatable :: analytic_gradient(:)
    real(dp), allocatable :: numeric_gradient(:)
    real(dp), allocatable :: gradient_error(:)
    real(dp), allocatable :: analytic_hessian(:, :)
    real(dp), allocatable :: numeric_hessian(:, :)
    real(dp), allocatable :: hessian_error(:, :)
    real(dp) :: max_gradient_error = huge(1.0_dp)
    real(dp) :: max_hessian_error = huge(1.0_dp)
    logical :: gradient_ok = .false.
    logical :: hessian_ok = .false.
    integer :: status = MAXLIK_INVALID_INPUT
  end type derivative_comparison

  public :: maxlik_objective, maxlik_gradient, maxlik_hessian, maxlik_scores

end module maxlik_types
