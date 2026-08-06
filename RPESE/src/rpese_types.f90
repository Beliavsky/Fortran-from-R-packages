! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_types
  use rpese_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: rpese_success = 0
  integer, parameter, public :: rpese_invalid_argument = 1
  integer, parameter, public :: rpese_numerical_failure = 2
  integer, parameter, public :: rpese_dependency_failure = 3
  integer, parameter, public :: rpese_unknown_estimator = 4

  integer, parameter, public :: se_if_iid = 1
  integer, parameter, public :: se_if_cor = 2
  integer, parameter, public :: se_if_cor_adapt = 3
  integer, parameter, public :: se_if_cor_pw = 4
  integer, parameter, public :: se_boot_iid = 5
  integer, parameter, public :: se_boot_cor = 6

  integer, parameter, public :: fit_exponential = 1
  integer, parameter, public :: fit_gamma = 2

  integer, parameter, public :: frequency_all = 1
  integer, parameter, public :: frequency_decimate = 2
  integer, parameter, public :: frequency_truncate = 3

  type, public :: rpese_options
    real(dp) :: alpha = 0.1_dp
    real(dp) :: beta = 0.1_dp
    real(dp) :: confidence = 0.95_dp
    real(dp) :: risk_free = 0.0_dp
    real(dp) :: threshold_constant = 0.0_dp
    integer :: lpm_order = 1
    character(len=8) :: sortino_threshold = 'mean'
    character(len=16) :: robust_family = 'mopt'
    real(dp) :: robust_efficiency = 0.95_dp
    logical :: clean_outliers = .false.
    integer :: fitting_method = fit_exponential
    integer :: polynomial_degree = 5
    real(dp) :: elastic_net_alpha = 0.5_dp
    integer :: frequency_mode = frequency_all
    real(dp) :: frequency_fraction = 0.5_dp
    real(dp) :: keep_fraction = 1.0_dp
    logical :: standardize_design = .false.
    real(dp) :: adaptive_a = 0.3_dp
    real(dp) :: adaptive_b = 0.7_dp
    integer :: bootstrap_replicates = 1000
    integer :: block_length = 0
    integer :: seed = 12345
    integer :: cv_folds = 5
    integer :: cv_repeats = 5
    integer :: num_lambda = 100
    integer :: max_iterations = 1000
    logical :: source_compatibility = .true.
  end type rpese_options

  type, public :: periodogram_result
    real(dp), allocatable :: frequency(:)
    real(dp), allocatable :: spectrum(:)
    integer :: status = rpese_success
    character(len=160) :: message = ''
  end type periodogram_result

  type, public :: se_result
    real(dp) :: estimate = 0.0_dp
    real(dp) :: standard_error = 0.0_dp
    real(dp) :: return_correlation = 0.0_dp
    real(dp) :: influence_correlation = 0.0_dp
    real(dp) :: prewhitened_influence_correlation = 0.0_dp
    real(dp) :: ar1_coefficient = 0.0_dp
    real(dp) :: adaptive_weight = 0.0_dp
    real(dp), allocatable :: coefficients(:)
    integer :: method = 0
    integer :: status = rpese_success
    character(len=32) :: estimator = ''
    character(len=160) :: message = ''
  end type se_result

  type, public :: matrix_se_result
    type(se_result), allocatable :: column(:)
    integer :: status = rpese_success
    character(len=160) :: message = ''
  end type matrix_se_result

end module rpese_types
