! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_types
  use robstattm_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: robstattm_success = 0
  integer, parameter, public :: robstattm_invalid_argument = 1
  integer, parameter, public :: robstattm_singular = 2
  integer, parameter, public :: robstattm_no_convergence = 3

  type, public :: robstattm_control
    character(len=16) :: family = 'bisquare'
    real(dp) :: bb = 0.5_dp
    real(dp) :: efficiency = 0.85_dp
    real(dp) :: tuning_chi = 1.54764_dp
    real(dp) :: tuning_psi = 3.4434_dp
    integer :: max_iter = 100
    integer :: n_resample = 500
    integer :: refine_iter = 50
    real(dp) :: tolerance = 1.0e-7_dp
    logical :: finite_sample_correction = .true.
  end type robstattm_control

  type, public :: location_scale_result
    real(dp) :: location = 0.0_dp
    real(dp) :: standard_error = 0.0_dp
    real(dp) :: scale = 0.0_dp
    integer :: iterations = 0
    integer :: status = robstattm_success
    logical :: converged = .false.
  end type location_scale_result

  type, public :: regression_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: standard_errors(:)
    real(dp) :: scale = 0.0_dp
    real(dp) :: objective = 0.0_dp
    real(dp) :: r_squared = 0.0_dp
    real(dp) :: adjusted_r_squared = 0.0_dp
    real(dp) :: mixing = 0.0_dp
    integer :: iterations = 0
    integer :: rank = 0
    integer :: status = robstattm_success
    logical :: converged = .false.
    character(len=24) :: method = ''
  end type regression_result

  type, public :: logistic_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residual_deviances(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: leverage_weights(:)
    real(dp) :: objective = 0.0_dp
    integer :: iterations = 0
    integer :: status = robstattm_success
    logical :: converged = .false.
    character(len=24) :: method = ''
  end type logistic_result

  type, public :: covariance_result
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: distances(:)
    real(dp), allocatable :: weights(:)
    integer, allocatable :: subset(:)
    real(dp) :: scale = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: objective = 0.0_dp
    integer :: iterations = 0
    integer :: status = robstattm_success
    logical :: converged = .false.
    character(len=48) :: method = ''
  end type covariance_result

  type, public :: pca_result
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: sdev(:)
    real(dp), allocatable :: loadings(:, :)
    real(dp), allocatable :: scores(:, :)
    real(dp), allocatable :: fitted(:, :)
    real(dp), allocatable :: initial_cumulative_variance(:)
    real(dp) :: explained_proportion = 0.0_dp
    integer :: n_components = 0
    integer :: iterations = 0
    integer :: status = robstattm_success
    logical :: converged = .false.
    character(len=32) :: method = ''
  end type pca_result

  type, public :: projection_result
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: distances(:)
    real(dp), allocatable :: outlyingness(:)
    logical, allocatable :: outliers(:)
    integer :: directions = 0
    integer :: status = robstattm_success
  end type projection_result

  type, public :: model_selection_result
    integer, allocatable :: selected_columns(:)
    real(dp), allocatable :: criterion_history(:)
    real(dp) :: criterion = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = robstattm_success
    logical :: converged = .false.
    character(len=24) :: direction = ''
  end type model_selection_result

  type, public :: linear_test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: chi_square_p_value = 1.0_dp
    real(dp) :: f_p_value = 1.0_dp
    integer :: df1 = 0
    integer :: df2 = 0
    integer :: status = robstattm_success
  end type linear_test_result
end module robstattm_types
