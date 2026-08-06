! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_types
  use ghyp_kinds, only : dp
  use tsgarch_types, only : garch_fit, garch_spec, fit_options
  implicit none
  private

  integer, parameter, public :: tsm_success = 0
  integer, parameter, public :: tsm_invalid_argument = 1
  integer, parameter, public :: tsm_no_convergence = 2
  integer, parameter, public :: tsm_numerical_failure = 3
  integer, parameter, public :: tsm_singular = 4

  type, public :: dcc_spec
    character(len=12) :: distribution = 'mvn'
    integer :: alpha_order = 1
    integer :: gamma_order = 0
    integer :: beta_order = 1
    logical :: constant_correlation = .false.
    integer :: correlation_method = 1
  end type dcc_spec

  type, public :: dcc_parameters
    real(dp), allocatable :: alpha(:)
    real(dp), allocatable :: gamma(:)
    real(dp), allocatable :: beta(:)
    real(dp) :: shape = 8.0_dp
  end type dcc_parameters

  type, public :: dcc_filter_result
    real(dp), allocatable :: q(:, :, :)
    real(dp), allocatable :: correlation(:, :, :)
    real(dp), allocatable :: covariance(:, :, :)
    real(dp), allocatable :: standardized_residuals(:, :)
    real(dp), allocatable :: sigma(:, :)
    real(dp), allocatable :: loglik_vector(:)
    real(dp), allocatable :: dcc_loglik_vector(:)
    real(dp), allocatable :: qbar(:, :)
    real(dp), allocatable :: nbar(:, :)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: dcc_log_likelihood = -huge(1.0_dp)
    integer :: status = tsm_invalid_argument
    character(len=240) :: message = ''
  end type dcc_filter_result

  type, public :: dcc_fit
    type(dcc_spec) :: spec
    type(dcc_parameters) :: parameters
    type(garch_fit), allocatable :: marginals(:)
    type(dcc_filter_result) :: filtered
    real(dp), allocatable :: packed_parameters(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: scores(:, :)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = tsm_invalid_argument
    character(len=240) :: message = ''
  end type dcc_fit

  type, public :: dcc_simulation
    real(dp), allocatable :: series(:, :, :)
    real(dp), allocatable :: sigma(:, :, :)
    real(dp), allocatable :: innovations(:, :, :)
    real(dp), allocatable :: correlation(:, :, :, :)
    integer :: status = tsm_invalid_argument
    character(len=200) :: message = ''
  end type dcc_simulation

  type, public :: dcc_forecast
    real(dp), allocatable :: mean(:, :)
    real(dp), allocatable :: covariance(:, :, :)
    real(dp), allocatable :: correlation(:, :, :)
    real(dp), allocatable :: sigma(:, :)
    real(dp), allocatable :: simulated(:, :, :)
    integer :: status = tsm_invalid_argument
    character(len=200) :: message = ''
  end type dcc_forecast

  type, public :: copula_spec
    character(len=12) :: distribution = 'gaussian'
    integer :: alpha_order = 1
    integer :: gamma_order = 0
    integer :: beta_order = 1
    logical :: constant_correlation = .false.
  end type copula_spec

  type, public :: copula_fit
    type(copula_spec) :: spec
    type(dcc_parameters) :: parameters
    type(garch_fit), allocatable :: marginals(:)
    type(dcc_filter_result) :: filtered
    real(dp), allocatable :: uniforms(:, :)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: status = tsm_invalid_argument
    character(len=240) :: message = ''
  end type copula_fit

  type, public :: ica_result
    real(dp), allocatable :: mixing(:, :)
    real(dp), allocatable :: unmixing(:, :)
    real(dp), allocatable :: whitening(:, :)
    real(dp), allocatable :: components(:, :)
    real(dp), allocatable :: mean(:)
    integer :: iterations = 0
    integer :: status = tsm_invalid_argument
    character(len=200) :: message = ''
  end type ica_result

  type, public :: gogarch_spec
    character(len=12) :: ica_method = 'fastica'
    type(garch_spec) :: factor_spec
  end type gogarch_spec

  type, public :: gogarch_fit
    type(gogarch_spec) :: spec
    type(ica_result) :: ica
    type(garch_fit), allocatable :: factors(:)
    real(dp), allocatable :: factor_series(:, :)
    real(dp), allocatable :: covariance(:, :, :)
    real(dp), allocatable :: correlation(:, :, :)
    real(dp), allocatable :: coskewness(:, :, :)
    real(dp), allocatable :: cokurtosis(:, :, :, :)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: status = tsm_invalid_argument
    character(len=240) :: message = ''
  end type gogarch_fit

  type, public :: gogarch_forecast
    real(dp), allocatable :: covariance(:, :, :)
    real(dp), allocatable :: correlation(:, :, :)
    real(dp), allocatable :: factor_variance(:, :)
    real(dp), allocatable :: simulated(:, :, :)
    integer :: status = tsm_invalid_argument
    character(len=200) :: message = ''
  end type gogarch_forecast

  type, public :: risk_result
    real(dp), allocatable :: probabilities(:)
    real(dp), allocatable :: value_at_risk(:, :)
    real(dp), allocatable :: expected_shortfall(:, :)
    integer :: status = tsm_invalid_argument
    character(len=160) :: message = ''
  end type risk_result


  type, public :: fft_distribution
    real(dp), allocatable :: grid(:)
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: cdf(:)
    integer :: status = tsm_invalid_argument
    character(len=160) :: message = ''
  end type fft_distribution

  type, public :: escc_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: degrees_freedom = 0
    integer :: status = tsm_invalid_argument
    character(len=160) :: message = ''
  end type escc_result

end module tsmarch_types
