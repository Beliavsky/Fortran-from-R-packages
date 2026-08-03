! SPDX-License-Identifier: GPL-3.0-only
module esback_types
  use esback_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: esback_ok = 0
  integer, parameter, public :: esback_invalid_input = 1
  integer, parameter, public :: esback_singular = 2
  integer, parameter, public :: esback_optimization_failed = 3
  integer, parameter, public :: esback_insufficient_data = 4
  integer, parameter, public :: esback_bootstrap_failed = 5

  integer, parameter, public :: sparsity_iid = 1
  integer, parameter, public :: sparsity_nid = 2
  integer, parameter, public :: sigma_ind = 1
  integer, parameter, public :: sigma_scl_n = 2
  integer, parameter, public :: sigma_scl_sp = 3
  integer, parameter, public :: bandwidth_bofinger = 1
  integer, parameter, public :: bandwidth_chamberlain = 2
  integer, parameter, public :: bandwidth_hall_sheather = 3

  type, public :: rng_state
    integer(kind=8) :: state = 1_8
  end type rng_state

  type, public :: er_backtest_result
    real(dp) :: pvalue_twosided_simple = 0.0_dp
    real(dp) :: pvalue_onesided_simple = 0.0_dp
    real(dp) :: pvalue_twosided_standardized = 0.0_dp
    real(dp) :: pvalue_onesided_standardized = 0.0_dp
    logical :: has_standardized = .false.
    integer :: n_exceedances = 0
    integer :: status = esback_ok
  end type er_backtest_result

  type, public :: cc_backtest_result
    real(dp) :: pvalue_twosided_simple = 0.0_dp
    real(dp) :: pvalue_onesided_simple = 0.0_dp
    real(dp) :: pvalue_twosided_general = 0.0_dp
    real(dp) :: pvalue_onesided_general = 0.0_dp
    logical :: has_general = .false.
    integer :: status = esback_ok
  end type cc_backtest_result

  type, public :: esreg_options
    integer :: sparsity = sparsity_nid
    integer :: sigma_est = sigma_scl_sp
    integer :: bandwidth_estimator = bandwidth_hall_sheather
    logical :: misspec = .true.
    integer :: early_stopping = 10
    integer :: max_iterations = 4000
    integer :: multistarts = 8
    real(dp) :: tolerance = 1.0e-8_dp
    integer(kind=8) :: seed = 104729_8
  end type esreg_options

  type, public :: esreg_fit_result
    real(dp), allocatable :: coefficients_q(:)
    real(dp), allocatable :: coefficients_e(:)
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: fitted_q(:)
    real(dp), allocatable :: fitted_e(:)
    real(dp) :: loss = 0.0_dp
    integer :: iterations = 0
    integer :: status = esback_ok
    logical :: covariance_available = .false.
  end type esreg_fit_result

  type, public :: esr_backtest_result
    real(dp) :: pvalue_twosided_asymptotic = 0.0_dp
    real(dp) :: pvalue_onesided_asymptotic = 0.0_dp
    real(dp) :: pvalue_twosided_bootstrap = 0.0_dp
    real(dp) :: pvalue_onesided_bootstrap = 0.0_dp
    real(dp) :: statistic = 0.0_dp
    integer :: version = 0
    integer :: bootstrap_successes = 0
    logical :: one_sided = .false.
    logical :: bootstrap_used = .false.
    type(esreg_fit_result) :: fit
    integer :: status = esback_ok
  end type esr_backtest_result
end module esback_types
