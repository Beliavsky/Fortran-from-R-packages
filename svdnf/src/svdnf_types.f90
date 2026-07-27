! SPDX-License-Identifier: GPL-3.0-only
module svdnf_types
  use svdnf_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: model_custom = 0
  integer, parameter, public :: model_heston = 1
  integer, parameter, public :: model_bates = 2
  integer, parameter, public :: model_duffie_pan_singleton = 3
  integer, parameter, public :: model_taylor = 4
  integer, parameter, public :: model_taylor_leverage = 5
  integer, parameter, public :: model_pitt_malik_doucet = 6
  integer, parameter, public :: model_capm_sv = 7

  abstract interface
    function state_function(x, parameters) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: parameters(:)
      real(dp) :: value
    end function state_function

    function jump_probability_function(n, parameters) result(value)
      import dp
      integer, intent(in) :: n
      real(dp), intent(in) :: parameters(:)
      real(dp) :: value
    end function jump_probability_function

    function jump_count_function(parameters) result(value)
      import dp
      real(dp), intent(in) :: parameters(:)
      integer :: value
    end function jump_count_function
  end interface

  type, public :: svm_dynamics
    integer :: model_id = model_heston
    character(len=32) :: model = 'Heston'
    real(dp) :: mu = 0.038_dp
    real(dp) :: kappa = 3.689_dp
    real(dp) :: theta = 0.032_dp
    real(dp) :: sigma = 0.446_dp
    real(dp) :: rho = -0.745_dp
    real(dp) :: omega = 5.125_dp
    real(dp) :: delta = 0.03_dp
    real(dp) :: alpha = -0.014_dp
    real(dp) :: rho_z = -1.809_dp
    real(dp) :: nu = 0.004_dp
    real(dp) :: p = 0.01_dp
    real(dp) :: phi = 0.965_dp
    real(dp) :: h = 1.0_dp / 252.0_dp
    real(dp), allocatable :: coefs(:)
    real(dp), allocatable :: mu_y_parameters(:)
    real(dp), allocatable :: sigma_y_parameters(:)
    real(dp), allocatable :: mu_x_parameters(:)
    real(dp), allocatable :: sigma_x_parameters(:)
    real(dp), allocatable :: jump_parameters(:)
    procedure(state_function), pointer, nopass :: custom_mu_y => null()
    procedure(state_function), pointer, nopass :: custom_sigma_y => null()
    procedure(state_function), pointer, nopass :: custom_mu_x => null()
    procedure(state_function), pointer, nopass :: custom_sigma_x => null()
    procedure(jump_probability_function), pointer, nopass :: custom_jump_probability => null()
    procedure(jump_count_function), pointer, nopass :: custom_jump_count => null()
  end type svm_dynamics

  type, public :: grid_type
    real(dp), allocatable :: var_mid_points(:)
    integer, allocatable :: jump_counts(:)
    real(dp), allocatable :: jump_mid_points(:)
  end type grid_type

  type, public :: simulation_result
    real(dp), allocatable :: volatility_factor(:)
    real(dp), allocatable :: returns(:)
    integer, allocatable :: jump_counts(:)
    real(dp), allocatable :: volatility_jumps(:)
    real(dp), allocatable :: return_jumps(:)
    logical :: ok = .true.
    character(len=160) :: message = ''
  end type simulation_result

  type, public :: filter_result
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp), allocatable :: likelihoods(:)
    real(dp), allocatable :: filter_grid(:,:)
    type(grid_type) :: grids
    type(svm_dynamics) :: dynamics
    real(dp), allocatable :: data(:)
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type filter_result

  type, public :: percentile_result
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: distributions(:,:)
    logical :: prediction = .false.
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type percentile_result

  type, public :: forecast_result
    real(dp), allocatable :: mean_volatility(:)
    real(dp), allocatable :: lower_volatility(:)
    real(dp), allocatable :: upper_volatility(:)
    real(dp), allocatable :: mean_return(:)
    real(dp), allocatable :: lower_return(:)
    real(dp), allocatable :: upper_return(:)
    real(dp) :: confidence = 0.95_dp
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type forecast_result

  type, public :: optimization_result
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
    logical :: ok = .false.
    character(len=160) :: message = ''
    type(filter_result) :: filter
  end type optimization_result

  public :: state_function, jump_probability_function, jump_count_function
end module svdnf_types
