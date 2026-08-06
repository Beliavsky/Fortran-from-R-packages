! SPDX-License-Identifier: GPL-3.0-only
module mass_types
  use rrcov_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: mass_success = 0
  integer, parameter, public :: mass_invalid_argument = 1
  integer, parameter, public :: mass_dimension_error = 2
  integer, parameter, public :: mass_singular = 3
  integer, parameter, public :: mass_no_convergence = 4
  integer, parameter, public :: mass_not_supported = 5

  type, public :: status_result
    integer :: status = mass_success
    integer :: iterations = 0
    character(len=160) :: message = ""
  end type status_result

  type, public :: huber_result
    real(dp) :: location = 0.0_dp
    real(dp) :: scale = 0.0_dp
    integer :: iterations = 0
    integer :: status = mass_success
  end type huber_result

  type, public :: regression_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: leverages(:)
    real(dp) :: sigma = 0.0_dp
    real(dp) :: log_likelihood = 0.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: theta = 0.0_dp
    integer :: rank = 0
    integer :: df_residual = 0
    integer :: iterations = 0
    integer :: status = mass_success
    character(len=48) :: method = ""
  end type regression_result

  type, public :: ridge_result
    real(dp), allocatable :: coefficients(:, :)
    real(dp), allocatable :: lambdas(:)
    real(dp), allocatable :: scales(:)
    real(dp), allocatable :: gcv(:)
    real(dp) :: y_mean = 0.0_dp
    integer :: status = mass_success
  end type ridge_result

  type, public :: density_fit_result
    real(dp), allocatable :: estimates(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp) :: log_likelihood = 0.0_dp
    real(dp) :: aic = 0.0_dp
    integer :: iterations = 0
    integer :: status = mass_success
    character(len=32) :: distribution = ""
  end type density_fit_result

  type, public :: ordinal_model
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: zeta(:)
    real(dp), allocatable :: covariance(:, :)
    integer :: n_categories = 0
    integer :: iterations = 0
    integer :: status = mass_success
    real(dp) :: log_likelihood = 0.0_dp
    real(dp) :: aic = 0.0_dp
    character(len=16) :: link = "logistic"
  end type ordinal_model

  type, public :: correspondence_result
    real(dp), allocatable :: correlations(:)
    real(dp), allocatable :: row_scores(:, :)
    real(dp), allocatable :: column_scores(:, :)
    integer :: status = mass_success
  end type correspondence_result

  type, public :: mca_result
    real(dp), allocatable :: singular_values(:)
    real(dp), allocatable :: row_scores(:, :)
    real(dp), allocatable :: factor_scores(:, :)
    integer :: status = mass_success
  end type mca_result

  type, public :: mds_result
    real(dp), allocatable :: points(:, :)
    real(dp), allocatable :: fitted_distances(:)
    real(dp) :: stress = 0.0_dp
    integer :: iterations = 0
    integer :: status = mass_success
    character(len=24) :: method = ""
  end type mds_result

  type, public :: kde2d_result
    real(dp), allocatable :: x_grid(:)
    real(dp), allocatable :: y_grid(:)
    real(dp), allocatable :: density(:, :)
    integer :: status = mass_success
  end type kde2d_result

  type, public :: loglinear_result
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: coefficients(:)
    real(dp) :: deviance = 0.0_dp
    real(dp) :: pearson = 0.0_dp
    integer :: df_residual = 0
    integer :: iterations = 0
    integer :: status = mass_success
  end type loglinear_result

  type, public :: model_selection_result
    logical, allocatable :: selected(:)
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: aic_path(:)
    integer, allocatable :: changed_column(:)
    integer :: steps = 0
    integer :: status = mass_success
  end type model_selection_result

end module mass_types
