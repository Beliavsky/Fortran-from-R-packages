! SPDX-License-Identifier: MIT
module bekks_types
  use bekks_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: bekk_full = 1
  integer, parameter, public :: bekk_diagonal = 2
  integer, parameter, public :: bekk_scalar = 3

  integer, parameter, public :: bekk_ok = 0
  integer, parameter, public :: bekk_invalid_input = 1
  integer, parameter, public :: bekk_invalid_parameters = 2
  integer, parameter, public :: bekk_linalg_failure = 3
  integer, parameter, public :: bekk_no_convergence = 4

  type, public :: bekk_spec_type
    integer :: model_type = bekk_full
    logical :: asymmetric = .false.
    real(dp), allocatable :: signs(:)
    real(dp), allocatable :: initial_theta(:)
  end type bekk_spec_type

  type, public :: bekk_parameters
    integer :: model_type = bekk_full
    logical :: asymmetric = .false.
    real(dp), allocatable :: c(:,:), a(:,:), b(:,:), g(:,:)
    real(dp) :: a_scalar = 0.0_dp
    real(dp) :: b_scalar = 0.0_dp
    real(dp) :: g_scalar = 0.0_dp
  end type bekk_parameters

  type, public :: bekk_fit_result
    type(bekk_spec_type) :: spec
    type(bekk_parameters) :: parameters
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: t_value(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: robust_covariance(:,:)
    real(dp), allocatable :: score(:,:)
    real(dp), allocatable :: h(:,:,:)
    real(dp), allocatable :: residuals(:,:)
    real(dp), allocatable :: likelihood_path(:)
    real(dp), allocatable :: data(:,:)
    real(dp), allocatable :: signs(:)
    real(dp) :: expected_indicator = 0.0_dp
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = bekk_ok
    logical :: stationary = .false.
    logical :: converged = .false.
  end type bekk_fit_result


  type, public :: bekk_filter_result
    real(dp), allocatable :: h(:,:,:)
    real(dp), allocatable :: residuals(:,:)
    integer :: status = bekk_ok
  end type bekk_filter_result

  type, public :: bekk_forecast_result
    real(dp), allocatable :: h(:,:,:)
    real(dp), allocatable :: covariance_lower(:,:,:)
    real(dp), allocatable :: covariance_upper(:,:,:)
    real(dp), allocatable :: standard_deviation(:,:)
    real(dp), allocatable :: correlation(:,:,:)
    integer :: status = bekk_ok
  end type bekk_forecast_result

  type, public :: bekk_var_result
    real(dp), allocatable :: value(:,:)
    real(dp), allocatable :: lower(:,:)
    real(dp), allocatable :: upper(:,:)
    real(dp), allocatable :: quantiles(:)
    integer :: status = bekk_ok
  end type bekk_var_result

  type, public :: bekk_backtest_result
    real(dp), allocatable :: var(:,:)
    real(dp), allocatable :: returns(:,:)
    real(dp), allocatable :: hit_rate(:)
    real(dp), allocatable :: kupiec_stat(:), kupiec_pvalue(:)
    real(dp), allocatable :: christoffersen_stat(:), christoffersen_pvalue(:)
    integer :: status = bekk_ok
  end type bekk_backtest_result

  type, public :: bekk_virf_result
    real(dp), allocatable :: response(:,:)
    real(dp), allocatable :: lower(:,:)
    real(dp), allocatable :: upper(:,:)
    integer :: status = bekk_ok
  end type bekk_virf_result

  type, public :: bekk_portmanteau_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: degrees_of_freedom = 0
    integer :: status = bekk_ok
  end type bekk_portmanteau_result

end module bekks_types
