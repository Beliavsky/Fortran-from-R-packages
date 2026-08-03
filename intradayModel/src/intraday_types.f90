! SPDX-License-Identifier: Apache-2.0
module intraday_types
  use intraday_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: intraday_ok = 0
  integer, parameter, public :: intraday_invalid_input = 1
  integer, parameter, public :: intraday_not_converged = 2
  integer, parameter, public :: intraday_numerical_failure = 3

  type, public :: volume_parameters
    real(dp) :: a_eta = 1.0_dp
    real(dp) :: a_mu = 0.0_dp
    real(dp) :: var_eta = 1.0e-4_dp
    real(dp) :: var_mu = 1.0e-4_dp
    real(dp) :: r = 1.0e-4_dp
    real(dp) :: x0(2) = [0.0_dp, 0.0_dp]
    real(dp) :: v0(2, 2) = reshape([1.0e-3_dp, 1.0e-7_dp, &
                                    1.0e-7_dp, 1.0e-5_dp], [2, 2])
    real(dp), allocatable :: phi(:)
  end type volume_parameters

  type, public :: parameter_mask
    logical :: a_eta = .false.
    logical :: a_mu = .false.
    logical :: var_eta = .false.
    logical :: var_mu = .false.
    logical :: r = .false.
    logical :: phi = .false.
    logical :: x0 = .false.
    logical :: v0 = .false.
  end type parameter_mask

  type, public :: volume_model_spec
    type(volume_parameters) :: fixed
    type(volume_parameters) :: initial
    type(parameter_mask) :: is_fixed
    type(parameter_mask) :: has_initial
  end type volume_model_spec

  type, public :: volume_fit_control
    logical :: acceleration = .true.
    integer :: maxit = 3000
    real(dp) :: abstol = 1.0e-4_dp
    logical :: save_history = .true.
    integer :: verbose = 0
  end type volume_fit_control

  type, public :: volume_model
    type(volume_parameters) :: par
    type(parameter_mask) :: fixed
    logical :: converged = .false.
    integer :: iterations = 0
    real(dp) :: final_change = huge(1.0_dp)
    integer :: status = intraday_ok
    character(len=160) :: message = ''
    type(volume_parameters), allocatable :: history(:)
  end type volume_model

  type, public :: error_metrics
    real(dp) :: mae = 0.0_dp
    real(dp) :: mape = 0.0_dp
    real(dp) :: rmse = 0.0_dp
  end type error_metrics

  type, public :: volume_decomposition
    real(dp), allocatable :: original_signal(:)
    real(dp), allocatable :: fitted_signal(:)
    real(dp), allocatable :: daily(:)
    real(dp), allocatable :: dynamic(:)
    real(dp), allocatable :: seasonal(:)
    real(dp), allocatable :: residual(:)
    type(error_metrics) :: error
    logical :: is_forecast = .false.
    integer :: status = intraday_ok
    character(len=160) :: message = ''
  end type volume_decomposition

  type, public :: kalman_output
    real(dp), allocatable :: x_pred(:, :)
    real(dp), allocatable :: v_pred(:, :, :)
    real(dp), allocatable :: x_filt(:, :)
    real(dp), allocatable :: v_filt(:, :, :)
    real(dp), allocatable :: gain(:, :)
    real(dp), allocatable :: x_smooth(:, :)
    real(dp), allocatable :: v_smooth(:, :, :)
    real(dp), allocatable :: smoother_gain(:, :, :)
    integer :: status = intraday_ok
    character(len=160) :: message = ''
  end type kalman_output

end module intraday_types
