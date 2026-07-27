! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_types
  use highfrequency_kinds, only: dp
  implicit none
  private

  type, public :: jump_test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp) :: critical_lower = 0.0_dp
    real(dp) :: critical_upper = 0.0_dp
    logical :: reject = .false.
  end type jump_test_result

  type, public :: iv_inference_result
    real(dp) :: estimate = 0.0_dp
    real(dp) :: standard_error = 0.0_dp
    real(dp) :: lower = 0.0_dp
    real(dp) :: upper = 0.0_dp
  end type iv_inference_result

  type, public :: har_model
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: sigma2 = 0.0_dp
    integer :: nobs = 0
    integer :: nreg = 0
    logical :: fitted_ok = .false.
  contains
    procedure :: predict_one => har_predict_one
  end type har_model

  type, public :: heavy_model
    real(dp) :: omega = 0.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: beta = 0.0_dp
    real(dp) :: omega_rm = 0.0_dp
    real(dp) :: alpha_rm = 0.0_dp
    real(dp) :: beta_rm = 0.0_dp
    real(dp), allocatable :: variance(:)
    real(dp), allocatable :: measure_variance(:)
    real(dp) :: loglik_returns = 0.0_dp
    real(dp) :: loglik_measure = 0.0_dp
    logical :: fitted_ok = .false.
  end type heavy_model

  type, public :: liquidity_result
    integer, allocatable :: direction(:)
    real(dp), allocatable :: midpoint(:)
    real(dp), allocatable :: effective_spread(:)
    real(dp), allocatable :: realized_spread(:)
    real(dp), allocatable :: price_impact(:)
    real(dp), allocatable :: proportional_effective_spread(:)
    real(dp), allocatable :: proportional_realized_spread(:)
    real(dp), allocatable :: proportional_price_impact(:)
    real(dp), allocatable :: depth_imbalance_difference(:)
    real(dp), allocatable :: depth_imbalance_ratio(:)
    real(dp), allocatable :: signed_trade_size(:)
    real(dp), allocatable :: quoted_spread(:)
    real(dp), allocatable :: proportional_quoted_spread(:)
  end type liquidity_result

  type, public :: lead_lag_result
    integer, allocatable :: lags(:)
    real(dp), allocatable :: contrast(:)
    integer :: optimal_lag = 0
    real(dp) :: maximum_contrast = 0.0_dp
  end type lead_lag_result

contains

  pure real(dp) function har_predict_one(self, regressors) result(value)
    class(har_model), intent(in) :: self
    real(dp), intent(in) :: regressors(:)
    if (.not. self%fitted_ok .or. size(regressors) + 1 /= size(self%coefficients)) then
      value = 0.0_dp
    else
      value = self%coefficients(1) + dot_product(self%coefficients(2:), regressors)
    end if
  end function har_predict_one

end module highfrequency_types
