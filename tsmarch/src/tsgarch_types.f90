! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_types
  use ghyp_kinds, only : dp
  use tsd_types, only : distribution_parameters
  implicit none
  private

  integer, parameter, public :: tsg_success = 0
  integer, parameter, public :: tsg_invalid_argument = 1
  integer, parameter, public :: tsg_no_convergence = 2
  integer, parameter, public :: tsg_numerical_failure = 3
  integer, parameter, public :: tsg_singular = 4

  integer, parameter, public :: init_unconditional = 1
  integer, parameter, public :: init_sample = 2
  integer, parameter, public :: init_backcast = 3

  type, public :: garch_spec
    character(len=12) :: model = 'garch'
    character(len=8) :: distribution = 'norm'
    integer :: p = 1
    integer :: q = 1
    logical :: constant = .false.
    logical :: variance_targeting = .false.
    logical :: multiplicative = .false.
    integer :: initialization = init_unconditional
    real(dp) :: backcast_lambda = 0.7_dp
    integer :: sample_n = 10
    real(dp) :: stationarity_limit = 0.999_dp
  end type garch_spec

  type, public :: garch_parameters
    real(dp) :: mu = 0.0_dp
    real(dp) :: omega = 0.01_dp
    real(dp) :: delta = 2.0_dp
    real(dp) :: rho = 0.90_dp
    real(dp) :: phi = 0.05_dp
    real(dp), allocatable :: alpha(:)
    real(dp), allocatable :: beta(:)
    real(dp), allocatable :: gamma(:)
    real(dp), allocatable :: eta(:)
    real(dp), allocatable :: xi(:)
    type(distribution_parameters) :: dist
  end type garch_parameters

  type, public :: fit_options
    integer :: max_iterations = 1800
    real(dp) :: tolerance = 1.0e-7_dp
    real(dp) :: simplex_scale = 0.08_dp
    logical :: compute_inference = .true.
  end type fit_options

  type, public :: garch_filter_result
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: variance(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: standardized_residuals(:)
    real(dp), allocatable :: loglik_vector(:)
    real(dp), allocatable :: permanent_component(:)
    real(dp), allocatable :: transitory_component(:)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: effective_omega = 0.0_dp
    real(dp) :: persistence = 0.0_dp
    real(dp) :: unconditional_variance = huge(1.0_dp)
    integer :: nobs = 0
    integer :: status = tsg_invalid_argument
    character(len=240) :: message = ''
  end type garch_filter_result

  type, public :: garch_fit
    type(garch_spec) :: spec
    type(garch_parameters) :: parameters
    type(garch_filter_result) :: filtered
    real(dp), allocatable :: packed_parameters(:)
    character(len=20), allocatable :: parameter_names(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: scores(:, :)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: npars = 0
    integer :: iterations = 0
    integer :: status = tsg_invalid_argument
    character(len=240) :: message = ''
  end type garch_fit

  type, public :: garch_forecast
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: variance(:)
    real(dp), allocatable :: probabilities(:)
    real(dp), allocatable :: quantiles(:, :)
    integer :: status = tsg_invalid_argument
    character(len=200) :: message = ''
  end type garch_forecast

  type, public :: garch_simulation
    real(dp), allocatable :: series(:, :)
    real(dp), allocatable :: sigma(:, :)
    real(dp), allocatable :: innovations(:, :)
    real(dp), allocatable :: permanent_component(:, :)
    real(dp), allocatable :: transitory_component(:, :)
    integer :: status = tsg_invalid_argument
    character(len=200) :: message = ''
  end type garch_simulation

  type, public :: backtest_result
    real(dp), allocatable :: actual(:)
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: value_at_risk(:)
    logical, allocatable :: exceedance(:)
    real(dp) :: expected_coverage = 0.0_dp
    real(dp) :: coverage = 0.0_dp
    real(dp) :: kupiec_statistic = 0.0_dp
    real(dp) :: kupiec_pvalue = 1.0_dp
    real(dp) :: independence_statistic = 0.0_dp
    real(dp) :: independence_pvalue = 1.0_dp
    real(dp) :: conditional_coverage_statistic = 0.0_dp
    real(dp) :: conditional_coverage_pvalue = 1.0_dp
    integer :: refits = 0
    integer :: status = tsg_invalid_argument
    character(len=200) :: message = ''
  end type backtest_result

  type, public :: profile_result
    character(len=20) :: parameter_name = ''
    real(dp), allocatable :: grid(:)
    real(dp), allocatable :: log_likelihood(:)
    integer :: status = tsg_invalid_argument
    character(len=200) :: message = ''
  end type profile_result

end module tsgarch_types
