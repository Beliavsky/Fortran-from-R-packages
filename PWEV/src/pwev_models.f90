! SPDX-License-Identifier: GPL-3.0-only
module pwev_models
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pwev_kinds, only : dp
  use pwev_status, only : PWEV_SUCCESS, PWEV_MODEL_FAILURE, PWEV_GARCH_MEAN, &
    PWEV_GARCH_SIGMA, PWEV_MEM_UPSTREAM_OOS, PWEV_MEM_RECURSIVE_OOS
  use pwev_types, only : pwev_control, PWEV_N_BASE_MODELS
  use rugarch_types, only : garch_fit_result
  use rugarch_fit, only : fit_garch11, fit_gjrgarch11, fit_igarch11
  use rugarch_models, only : forecast_volatility
  use rumidas_types, only : mem_spec, rumidas_fit_control, rumidas_fit_result, RUMIDAS_MEM
  use rumidas_status, only : RUMIDAS_SUCCESS
  use rumidas_fit, only : umemfit
  use rumidas_mem_models, only : mem_pred_no_skew
  implicit none
  private
  public :: fit_pwev_base_models
contains

  subroutine fit_pwev_base_models(train_data, test_data, control, train_models, test_models, model_status)
    real(dp), intent(in) :: train_data(:), test_data(:)
    type(pwev_control), intent(in) :: control
    real(dp), intent(out) :: train_models(size(train_data), PWEV_N_BASE_MODELS)
    real(dp), intent(out) :: test_models(size(test_data), PWEV_N_BASE_MODELS)
    integer, intent(out) :: model_status(PWEV_N_BASE_MODELS)
    type(garch_fit_result) :: fit

    fit = fit_garch11(train_data, fit_mean=.true., max_iterations=control%garch_max_iterations)
    call extract_garch_predictions(fit, train_data, size(test_data), control%garch_output, &
      train_models(:, 1), test_models(:, 1), model_status(1))

    fit = fit_gjrgarch11(train_data, fit_mean=.true., max_iterations=control%garch_max_iterations)
    call extract_garch_predictions(fit, train_data, size(test_data), control%garch_output, &
      train_models(:, 2), test_models(:, 2), model_status(2))

    fit = fit_igarch11(train_data, fit_mean=.true., max_iterations=control%garch_max_iterations)
    call extract_garch_predictions(fit, train_data, size(test_data), control%garch_output, &
      train_models(:, 3), test_models(:, 3), model_status(3))

    call fit_mem_predictions(train_data, test_data, control, train_models(:, 4), test_models(:, 4), &
      model_status(4))
  end subroutine fit_pwev_base_models

  subroutine extract_garch_predictions(fit, train_data, horizon, output_mode, train_prediction, &
      test_prediction, status)
    type(garch_fit_result), intent(in) :: fit
    real(dp), intent(in) :: train_data(:)
    integer, intent(in) :: horizon, output_mode
    real(dp), intent(out) :: train_prediction(size(train_data)), test_prediction(horizon)
    integer, intent(out) :: status
    real(dp) :: fallback

    fallback = sum(train_data) / real(size(train_data), dp)
    train_prediction = fallback
    test_prediction = fallback
    status = PWEV_MODEL_FAILURE
    if (.not. allocated(fit%residuals) .or. .not. allocated(fit%sigma)) return
    if (size(fit%residuals) /= size(train_data) .or. size(fit%sigma) /= size(train_data)) return
    if (any(.not. ieee_is_finite(fit%residuals)) .or. any(.not. ieee_is_finite(fit%sigma))) return

    select case (output_mode)
    case (PWEV_GARCH_MEAN)
      train_prediction = train_data - fit%residuals
      test_prediction = fit%spec%mean
    case (PWEV_GARCH_SIGMA)
      train_prediction = fit%sigma
      call forecast_volatility(fit%spec, fit%residuals, fit%sigma, horizon, test_prediction)
    case default
      return
    end select
    if (all(ieee_is_finite(train_prediction)) .and. all(ieee_is_finite(test_prediction))) &
      status = PWEV_SUCCESS
  end subroutine extract_garch_predictions

  subroutine fit_mem_predictions(train_data, test_data, control, train_prediction, test_prediction, status)
    real(dp), intent(in) :: train_data(:), test_data(:)
    type(pwev_control), intent(in) :: control
    real(dp), intent(out) :: train_prediction(size(train_data)), test_prediction(size(test_data))
    integer, intent(out) :: status
    type(mem_spec) :: spec
    type(rumidas_fit_control) :: fit_control
    type(rumidas_fit_result) :: fit
    real(dp), allocatable :: upstream_prediction(:)
    real(dp) :: fallback
    integer :: rumidas_status

    fallback = sum(train_data) / real(size(train_data), dp)
    train_prediction = fallback
    test_prediction = fallback
    status = PWEV_MODEL_FAILURE
    if (any(train_data <= 0.0_dp) .or. any(test_data <= 0.0_dp)) return

    spec%model = RUMIDAS_MEM
    spec%skew = .false.
    spec%k = 0
    fit_control%random_starts = max(1, control%mem_random_starts)
    fit_control%max_iterations = max(1, control%mem_max_iterations)
    fit_control%random_seed = control%random_seed
    fit_control%method = 'bfgs'
    fit_control%compute_robust_covariance = .false.
    call umemfit(spec, train_data, fit, rumidas_status, control=fit_control)
    if (rumidas_status /= RUMIDAS_SUCCESS .or. .not. fit%converged) return
    if (.not. allocated(fit%conditional) .or. .not. allocated(fit%coefficients)) return
    if (size(fit%conditional) /= size(train_data) .or. size(fit%coefficients) < 2) return
    train_prediction = fit%conditional

    select case (control%mem_oos_mode)
    case (PWEV_MEM_UPSTREAM_OOS)
      call mem_pred_no_skew(fit%coefficients, test_data, upstream_prediction, rumidas_status)
      if (rumidas_status /= RUMIDAS_SUCCESS .or. size(upstream_prediction) /= size(test_data)) return
      test_prediction = upstream_prediction
    case (PWEV_MEM_RECURSIVE_OOS)
      call recursive_mem_forecast(train_data, fit%conditional, fit%coefficients, test_prediction)
    case default
      return
    end select
    if (all(ieee_is_finite(train_prediction)) .and. all(ieee_is_finite(test_prediction)) .and. &
        all(train_prediction > 0.0_dp) .and. all(test_prediction > 0.0_dp)) status = PWEV_SUCCESS
  end subroutine fit_mem_predictions

  subroutine recursive_mem_forecast(train_data, fitted, coefficients, forecast)
    real(dp), intent(in) :: train_data(:), fitted(:), coefficients(:)
    real(dp), intent(out) :: forecast(:)
    real(dp) :: alpha, beta, intercept, mean_train, last_x, last_short
    integer :: h

    alpha = coefficients(1)
    beta = coefficients(2)
    mean_train = sum(train_data) / real(size(train_data), dp)
    intercept = (1.0_dp - alpha - beta) * mean_train
    last_x = train_data(size(train_data))
    last_short = fitted(size(fitted))
    do h = 1, size(forecast)
      forecast(h) = intercept + alpha * last_x + beta * last_short
      forecast(h) = max(forecast(h), tiny(1.0_dp))
      last_x = forecast(h)
      last_short = forecast(h)
    end do
  end subroutine recursive_mem_forecast

end module pwev_models
