! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_types
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error
  implicit none
  private

  type, public :: local_pca_point_result
    real(dp), allocatable :: factors(:,:)
    real(dp), allocatable :: f_hat(:)
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: weights(:)
    type(tvmvp_error) :: error
  end type local_pca_point_result

  type, public :: local_pca_result
    real(dp), allocatable :: factors(:,:,:)
    real(dp), allocatable :: loadings(:,:,:)
    real(dp), allocatable :: weights(:,:)
    real(dp), allocatable :: f_hat(:,:)
    integer :: m = 0
    type(tvmvp_error) :: error
  end type local_pca_result

  type, public :: factor_selection_result
    integer :: optimal_m = 0
    real(dp), allocatable :: ic_values(:)
    real(dp), allocatable :: residual_variances(:)
    real(dp), allocatable :: penalties(:)
    type(tvmvp_error) :: error
  end type factor_selection_result

  type, public :: poet_result
    real(dp) :: best_rho = 0.0_dp
    real(dp) :: rho_lower = 0.0_dp
    real(dp) :: min_frobenius = 0.0_dp
    real(dp), allocatable :: residual_cov(:,:)
    real(dp), allocatable :: total_cov(:,:)
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: naive_residual_cov(:,:)
    type(tvmvp_error) :: error
  end type poet_result

  type, public :: hypothesis_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp) :: m_hat = 0.0_dp
    real(dp) :: bias_term = 0.0_dp
    real(dp) :: variance_term = 0.0_dp
    real(dp), allocatable :: bootstrap_statistics(:)
    type(tvmvp_error) :: error
  end type hypothesis_result

  type, public :: portfolio_result
    real(dp), allocatable :: weights(:)
    real(dp) :: expected_return = 0.0_dp
    real(dp) :: risk = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    logical :: available = .false.
    type(tvmvp_error) :: error
  end type portfolio_result

  type, public :: portfolio_prediction_result
    type(portfolio_result) :: minimum_variance
    type(portfolio_result) :: maximum_sharpe
    type(portfolio_result) :: return_constrained
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: expected_asset_returns(:)
    integer :: factors = 0
    real(dp) :: bandwidth = 0.0_dp
    type(tvmvp_error) :: error
  end type portfolio_prediction_result

  type, public :: performance_metrics
    real(dp) :: cumulative_excess_return = 0.0_dp
    real(dp) :: mean_excess_return = 0.0_dp
    real(dp) :: standard_deviation = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: annualized_mean = 0.0_dp
    real(dp) :: annualized_standard_deviation = 0.0_dp
  end type performance_metrics

  type, public :: expanding_window_result
    integer, allocatable :: rebalance_dates(:)
    integer, allocatable :: holding_lengths(:)
    real(dp), allocatable :: weights(:,:)
    real(dp), allocatable :: tvmvp_returns(:)
    real(dp), allocatable :: equal_returns(:)
    type(performance_metrics) :: tvmvp_metrics
    type(performance_metrics) :: equal_metrics
    type(tvmvp_error) :: error
  end type expanding_window_result
end module tvmvp_types
