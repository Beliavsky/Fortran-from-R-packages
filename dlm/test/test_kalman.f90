! SPDX-License-Identifier: GPL-2.0-or-later
program test_kalman
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use dlm, only : dp, dlm_model, dlm_filter_result, dlm_smooth_result, dlm_forecast_result
    use dlm, only : dlm_success, dlm_filter, dlm_forecast, dlm_ll, dlm_mod_poly, dlm_residuals, dlm_smooth
    implicit none

    type(dlm_model) :: model, future
    type(dlm_filter_result) :: filtered
    type(dlm_smooth_result) :: smoothed
    type(dlm_forecast_result) :: forecast
    real(dp) :: y(5, 1), nll
    real(dp), allocatable :: residual(:, :), residual_sd(:, :)
    integer :: info

    call dlm_mod_poly(1, model, info, d_v=0.5_dp, d_w=[0.1_dp], m0=[0.0_dp], &
                      c0=reshape([1.0_dp], [1, 1]))
    if (info /= dlm_success) error stop "local-level construction failed"
    y(:, 1) = [1.0_dp, 2.0_dp, 1.5_dp, 0.0_dp, 1.2_dp]
    y(4, 1) = ieee_value(0.0_dp, ieee_quiet_nan)

    call dlm_filter(y, model, filtered, info)
    if (info /= dlm_success) error stop "dlmFilter failed"
    if (abs(filtered%m(2, 1) - 0.6875_dp) > 2.0e-12_dp) error stop "first filtered mean mismatch"
    if (abs(filtered%m(6, 1) - 1.3016113410231576_dp) > 2.0e-12_dp) error stop "final filtered mean mismatch"
    if (abs(filtered%c(1, 1, 6) - 0.22241789204895662_dp) > 2.0e-12_dp) error stop "final filtered variance mismatch"
    if (abs(filtered%m(5, 1) - filtered%a(4, 1)) > 2.0e-12_dp) error stop "all-missing update mismatch"

    call dlm_ll(y, model, nll, info)
    if (info /= dlm_success) error stop "dlmLL failed"
    if (abs(nll - 1.3302399037948809_dp) > 2.0e-11_dp) error stop "dlmLL value mismatch"

    call dlm_smooth(filtered, smoothed, info)
    if (info /= dlm_success) error stop "dlmSmooth failed"
    if (abs(smoothed%s(1, 1) - 1.0781896551724138_dp) > 2.0e-8_dp) error stop "smoothed initial mean mismatch"
    if (abs(smoothed%s(6, 1) - filtered%m(6, 1)) > 2.0e-12_dp) error stop "smoothed terminal state mismatch"

    call dlm_residuals(filtered, residual, residual_sd)
    if (.not. ieee_is_nan(residual(4, 1))) error stop "missing residual was not preserved"
    if (abs(residual_sd(1, 1) - sqrt(1.6_dp)) > 2.0e-12_dp) error stop "forecast residual SD mismatch"

    future = model
    future%m0 = filtered%m(6, :)
    future%c0 = filtered%c(:, :, 6)
    call dlm_forecast(future, 2, forecast, info)
    if (info /= dlm_success) error stop "dlmForecast failed"
    if (abs(forecast%a(1, 1) - filtered%m(6, 1)) > 2.0e-12_dp) error stop "forecast mean mismatch"
    if (abs(forecast%q(1, 1, 1) - 0.8224178920489567_dp) > 2.0e-12_dp) error stop "forecast variance mismatch"

    print *, "test_kalman: PASS"
end program test_kalman
