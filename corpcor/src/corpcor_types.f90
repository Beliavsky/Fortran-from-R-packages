! SPDX-License-Identifier: GPL-3.0-or-later
module corpcor_types
  use corpcor_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: corpcor_success = 0
  integer, parameter, public :: corpcor_invalid_argument = 1
  integer, parameter, public :: corpcor_dimension_error = 2
  integer, parameter, public :: corpcor_numerical_error = 3

  type, public :: svd_result
    real(dp), allocatable :: d(:)
    real(dp), allocatable :: u(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp) :: tol = 0.0_dp
    integer :: rank = 0
    integer :: status = corpcor_success
  end type svd_result

  type, public :: rank_condition_result
    integer :: rank = 0
    real(dp) :: condition = huge(1.0_dp)
    real(dp) :: tol = 0.0_dp
    integer :: status = corpcor_success
  end type rank_condition_result

  type, public :: moments_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: variance(:)
    integer :: status = corpcor_success
  end type moments_result

  type, public :: scale_result
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
    logical, allocatable :: zero_scale(:)
    integer :: status = corpcor_success
  end type scale_result

  type, public :: matrix_shrinkage_result
    real(dp), allocatable :: value(:, :)
    real(dp) :: lambda = 0.0_dp
    logical :: lambda_estimated = .false.
    real(dp) :: lambda_var = 0.0_dp
    logical :: lambda_var_estimated = .false.
    real(dp), allocatable :: standardized_partial_variance(:)
    integer :: status = corpcor_success
  end type matrix_shrinkage_result

  type, public :: vector_shrinkage_result
    real(dp), allocatable :: value(:)
    real(dp) :: lambda_var = 0.0_dp
    logical :: lambda_var_estimated = .false.
    real(dp), allocatable :: standardized_partial_variance(:)
    integer :: status = corpcor_success
  end type vector_shrinkage_result

  type, public :: covariance_decomposition
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: variance(:)
    integer :: status = corpcor_success
  end type covariance_decomposition

  type, public :: precision_decomposition
    real(dp), allocatable :: partial_correlation(:, :)
    real(dp), allocatable :: partial_variance(:)
    integer :: status = corpcor_success
  end type precision_decomposition
end module corpcor_types
