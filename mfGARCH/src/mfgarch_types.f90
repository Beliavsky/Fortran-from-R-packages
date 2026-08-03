! SPDX-License-Identifier: MIT
module mfgarch_types
  use mfgarch_kinds, only : dp
  implicit none
  private

  type, public :: mfgarch_model
    real(dp) :: mu = 0.0_dp
    real(dp) :: alpha = 0.02_dp
    real(dp) :: beta = 0.85_dp
    real(dp) :: gamma = 0.04_dp
    real(dp) :: m = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: w1 = 1.0_dp
    real(dp) :: w2 = 3.0_dp
    real(dp) :: theta_two = 0.0_dp
    real(dp) :: w1_two = 1.0_dp
    real(dp) :: w2_two = 3.0_dp
    integer :: k = 0
    integer :: k_two = 0
    logical :: asymmetric = .true.
    logical :: unrestricted_weights = .false.
    logical :: has_second = .false.
  end type mfgarch_model

  type, public :: mfgarch_fit_control
    integer :: max_iterations = 500
    integer :: max_function_evaluations = 20000
    real(dp) :: tolerance = 1.0e-7_dp
    real(dp) :: gradient_step = 1.0e-5_dp
    real(dp) :: initial_simplex_step = 0.15_dp
    logical :: multi_start = .true.
    logical :: compute_inference = .true.
    logical :: trace = .false.
    character(len=16) :: method = 'bfgs'
  end type mfgarch_fit_control

  type, public :: mfgarch_fit_result
    type(mfgarch_model) :: model
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    real(dp) :: variance_ratio = 0.0_dp
    real(dp) :: tau_forecast = 0.0_dp
    integer :: status = 0
    integer :: iterations = 0
    integer :: function_evaluations = 0
    character(len=160) :: message = ''
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: g(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: weights_two(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: robust_covariance(:,:)
    real(dp), allocatable :: opg_covariance(:,:)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: robust_standard_error(:)
    real(dp), allocatable :: opg_standard_error(:)
  end type mfgarch_fit_result

  type, public :: mfgarch_simulation
    real(dp), allocatable :: returns(:)
    real(dp), allocatable :: covariate(:)
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: g(:)
    real(dp), allocatable :: realized_variance(:)
    real(dp), allocatable :: realized_variance_half_hour(:)
    real(dp), allocatable :: rv_5(:)
    real(dp), allocatable :: rv_22(:)
    real(dp), allocatable :: rv_half_hour_5(:)
    real(dp), allocatable :: rv_half_hour_22(:)
    real(dp), allocatable :: intraday_returns(:)
    integer, allocatable :: low_frequency_period(:)
  end type mfgarch_simulation

end module mfgarch_types
