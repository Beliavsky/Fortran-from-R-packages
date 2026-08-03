! SPDX-License-Identifier: GPL-3.0-or-later
module corpcor
  use corpcor_kinds, only : dp
  use corpcor_types
  use corpcor_linalg, only : fast_svd, pseudoinverse, matrix_power, rank_condition, &
    is_positive_definite, make_positive_definite, covariance_to_correlation
  use corpcor_weighted, only : normalized_weights, weighted_variance, weighted_moments, &
    weighted_scale, median_value
  use corpcor_matrix_tools, only : symmetric_matrix_to_vector, symmetric_matrix_indices, &
    vector_to_symmetric_matrix, rebuild_covariance, decompose_covariance, &
    rebuild_precision, decompose_precision
  use corpcor_shrinkage, only : estimate_lambda, estimate_lambda_var, variance_shrinkage, &
    correlation_power_shrinkage, correlation_shrinkage, inverse_correlation_shrinkage, &
    covariance_shrinkage, inverse_covariance_shrinkage, partial_correlation_shrinkage, &
    partial_variance_shrinkage, correlation_to_partial, partial_to_correlation, &
    crossprod_correlation_power_shrinkage
  implicit none
  private

  public :: dp
  public :: svd_result, rank_condition_result, moments_result, scale_result
  public :: matrix_shrinkage_result, vector_shrinkage_result
  public :: covariance_decomposition, precision_decomposition
  public :: corpcor_success, corpcor_invalid_argument, corpcor_dimension_error
  public :: corpcor_numerical_error
  public :: fast_svd, pseudoinverse, matrix_power, rank_condition
  public :: is_positive_definite, make_positive_definite, covariance_to_correlation
  public :: normalized_weights, weighted_variance, weighted_moments, weighted_scale
  public :: median_value
  public :: symmetric_matrix_to_vector, symmetric_matrix_indices
  public :: vector_to_symmetric_matrix, rebuild_covariance, decompose_covariance
  public :: rebuild_precision, decompose_precision
  public :: estimate_lambda, estimate_lambda_var, variance_shrinkage
  public :: correlation_power_shrinkage, correlation_shrinkage
  public :: inverse_correlation_shrinkage, covariance_shrinkage
  public :: inverse_covariance_shrinkage, partial_correlation_shrinkage
  public :: partial_variance_shrinkage, correlation_to_partial
  public :: partial_to_correlation, crossprod_correlation_power_shrinkage

  ! Compact aliases corresponding to the upstream R names with dots replaced by underscores.
  public :: wt_var, wt_moments, wt_scale, mpower
  public :: sm2vec, sm_index, vec2sm, rebuild_cov, decompose_cov
  public :: rebuild_invcov, decompose_invcov, cor2pcor, pcor2cor
  public :: var_shrink, cor_shrink, invcor_shrink, cov_shrink, invcov_shrink
  public :: pcor_shrink, pvar_shrink, powcor_shrink, crossprod_powcor_shrink

contains

  function wt_var(x, w, status) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: w(:)
    integer, intent(out), optional :: status
    real(dp) :: v
    v = weighted_variance(x, w, status)
  end function wt_var

  function wt_moments(x, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    type(moments_result) :: res
    res = weighted_moments(x, w)
  end function wt_moments

  function wt_scale(x, w, center, scale) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: w(:)
    logical, intent(in), optional :: center, scale
    type(scale_result) :: res
    res = weighted_scale(x, w, center, scale)
  end function wt_scale

  function mpower(a, alpha, pseudo, tolerance, status) result(out)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: alpha
    logical, intent(in), optional :: pseudo
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: out(:, :)
    out = matrix_power(a, alpha, pseudo, tolerance, status)
  end function mpower

  function sm2vec(a, diag, status) result(v)
    real(dp), intent(in) :: a(:, :)
    logical, intent(in), optional :: diag
    integer, intent(out), optional :: status
    real(dp), allocatable :: v(:)
    v = symmetric_matrix_to_vector(a, diag, status)
  end function sm2vec

  function sm_index(n, diag, status) result(index)
    integer, intent(in) :: n
    logical, intent(in), optional :: diag
    integer, intent(out), optional :: status
    integer, allocatable :: index(:, :)
    index = symmetric_matrix_indices(n, diag, status)
  end function sm_index

  function vec2sm(v, diag, order, status) result(a)
    real(dp), intent(in) :: v(:)
    logical, intent(in), optional :: diag
    integer, intent(in), optional :: order(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: a(:, :)
    a = vector_to_symmetric_matrix(v, diag, order, status)
  end function vec2sm

  function rebuild_cov(r, v, status) result(a)
    real(dp), intent(in) :: r(:, :), v(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: a(:, :)
    a = rebuild_covariance(r, v, status)
  end function rebuild_cov

  function decompose_cov(a) result(res)
    real(dp), intent(in) :: a(:, :)
    type(covariance_decomposition) :: res
    res = decompose_covariance(a)
  end function decompose_cov

  function rebuild_invcov(pr, pv, status) result(a)
    real(dp), intent(in) :: pr(:, :), pv(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: a(:, :)
    a = rebuild_precision(pr, pv, status)
  end function rebuild_invcov

  function decompose_invcov(a) result(res)
    real(dp), intent(in) :: a(:, :)
    type(precision_decomposition) :: res
    res = decompose_precision(a)
  end function decompose_invcov

  function cor2pcor(a, tolerance, status) result(out)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: out(:, :)
    out = correlation_to_partial(a, tolerance, status)
  end function cor2pcor

  function pcor2cor(a, tolerance, status) result(out)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: out(:, :)
    out = partial_to_correlation(a, tolerance, status)
  end function pcor2cor

  function var_shrink(x, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda_var, w(:)
    type(vector_shrinkage_result) :: res
    res = variance_shrinkage(x, lambda_var, w)
  end function var_shrink

  function cor_shrink(x, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, w(:)
    type(matrix_shrinkage_result) :: res
    res = correlation_shrinkage(x, lambda, w)
  end function cor_shrink

  function invcor_shrink(x, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, w(:)
    type(matrix_shrinkage_result) :: res
    res = inverse_correlation_shrinkage(x, lambda, w)
  end function invcor_shrink

  function cov_shrink(x, lambda, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, lambda_var, w(:)
    type(matrix_shrinkage_result) :: res
    res = covariance_shrinkage(x, lambda, lambda_var, w)
  end function cov_shrink

  function invcov_shrink(x, lambda, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, lambda_var, w(:)
    type(matrix_shrinkage_result) :: res
    res = inverse_covariance_shrinkage(x, lambda, lambda_var, w)
  end function invcov_shrink

  function pcor_shrink(x, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, w(:)
    type(matrix_shrinkage_result) :: res
    res = partial_correlation_shrinkage(x, lambda, w)
  end function pcor_shrink

  function pvar_shrink(x, lambda, lambda_var, w) result(res)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: lambda, lambda_var, w(:)
    type(vector_shrinkage_result) :: res
    res = partial_variance_shrinkage(x, lambda, lambda_var, w)
  end function pvar_shrink

  function powcor_shrink(x, alpha, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :), alpha
    real(dp), intent(in), optional :: lambda, w(:)
    type(matrix_shrinkage_result) :: res
    res = correlation_power_shrinkage(x, alpha, lambda, w)
  end function powcor_shrink

  function crossprod_powcor_shrink(x, y, alpha, lambda, w) result(res)
    real(dp), intent(in) :: x(:, :), y(:, :), alpha
    real(dp), intent(in), optional :: lambda, w(:)
    type(matrix_shrinkage_result) :: res
    res = crossprod_correlation_power_shrinkage(x, y, alpha, lambda, w)
  end function crossprod_powcor_shrink

end module corpcor
