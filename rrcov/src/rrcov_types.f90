! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_types
  use rrcov_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: rrcov_success = 0
  integer, parameter, public :: rrcov_invalid_argument = 1
  integer, parameter, public :: rrcov_dimension_error = 2
  integer, parameter, public :: rrcov_singular = 3
  integer, parameter, public :: rrcov_no_convergence = 4
  integer, parameter, public :: rrcov_allocation_error = 5

  type, public :: covariance_result
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: distances(:)
    real(dp), allocatable :: weights(:)
    integer, allocatable :: subset(:)
    integer :: n_obs = 0
    integer :: rank = 0
    integer :: iterations = 0
    integer :: status = rrcov_success
    real(dp) :: objective = 0.0_dp
    character(len=48) :: method = ""
  end type covariance_result

  type, public :: pca_result
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: loadings(:, :)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: scores(:, :)
    real(dp), allocatable :: score_distances(:)
    real(dp), allocatable :: orthogonal_distances(:)
    integer :: n_obs = 0
    integer :: n_components = 0
    integer :: status = rrcov_success
    character(len=48) :: method = ""
  end type pca_result

  type, public :: lda_model
    real(dp), allocatable :: means(:, :)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: inverse(:, :)
    real(dp), allocatable :: priors(:)
    real(dp), allocatable :: coefficients(:, :)
    real(dp), allocatable :: constants(:)
    integer, allocatable :: labels(:)
    integer :: n_features = 0
    integer :: n_groups = 0
    integer :: status = rrcov_success
    character(len=48) :: method = ""
  end type lda_model

  type, public :: qda_model
    real(dp), allocatable :: means(:, :)
    real(dp), allocatable :: covariance(:, :, :)
    real(dp), allocatable :: inverse(:, :, :)
    real(dp), allocatable :: log_determinant(:)
    real(dp), allocatable :: priors(:)
    integer, allocatable :: labels(:)
    integer :: n_features = 0
    integer :: n_groups = 0
    integer :: status = rrcov_success
    character(len=48) :: method = ""
  end type qda_model

  type, public :: test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp) :: df1 = 0.0_dp
    real(dp) :: df2 = 0.0_dp
    real(dp) :: lambda = 1.0_dp
    integer :: status = rrcov_success
    character(len=64) :: method = ""
  end type test_result

end module rrcov_types
