! SPDX-License-Identifier: GPL-3.0-only
module rumidas_forecast
  use rumidas_kinds, only: dp
  use rumidas_status
  use rumidas_types
  implicit none
  private
  public :: multi_step_ahead_pred, garch_midas_multi_step_forecast

contains

  subroutine garch_midas_multi_step_forecast(param, spec, last_short_run_variance, &
      last_long_run_variance, horizon, forecast, status, x_last, x_ar1)
    real(dp), intent(in) :: param(:)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: last_short_run_variance, last_long_run_variance
    integer, intent(in) :: horizon
    real(dp), allocatable, intent(out) :: forecast(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: x_last, x_ar1

    real(dp) :: alpha, beta, gamma_value, persistence, x_coefficient, delta, x_effect
    integer :: i, j, index
    logical :: has_x

    status = RUMIDAS_SUCCESS
    allocate(forecast(max(horizon, 0)))
    if (horizon <= 0 .or. last_short_run_variance <= 0.0_dp .or. last_long_run_variance <= 0.0_dp) then
      status = RUMIDAS_INVALID_INPUT
      return
    end if
    alpha = param(1)
    beta = param(2)
    index = 3
    gamma_value = 0.0_dp
    if (spec%skew) then
      gamma_value = param(index)
      index = index + 1
    end if
    has_x = spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX
    x_coefficient = 0.0_dp
    if (has_x) then
      if (.not. present(x_last) .or. .not. present(x_ar1)) then
        status = RUMIDAS_INVALID_INPUT
        forecast = 0.0_dp
        return
      end if
      x_coefficient = param(index)
    end if
    persistence = alpha + beta + 0.5_dp * gamma_value
    if (persistence >= 1.0_dp .or. persistence < 0.0_dp) then
      status = RUMIDAS_INVALID_PARAMETER
      forecast = 0.0_dp
      return
    end if

    forecast(1) = last_short_run_variance * last_long_run_variance
    do i = 2, horizon
      x_effect = 0.0_dp
      if (has_x) then
        delta = x_ar1
        do j = 0, i - 2
          x_effect = x_effect + persistence ** j * delta ** (i - j - 1) * x_coefficient * x_last
        end do
      end if
      forecast(i) = (1.0_dp + persistence ** (i - 1) * (last_short_run_variance - 1.0_dp) + x_effect) * &
        last_long_run_variance
    end do
  end subroutine garch_midas_multi_step_forecast

  subroutine multi_step_ahead_pred(param, spec, last_short_run_variance, last_long_run_variance, &
      horizon, forecast, status, x_last, x_ar1)
    real(dp), intent(in) :: param(:)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: last_short_run_variance, last_long_run_variance
    integer, intent(in) :: horizon
    real(dp), allocatable, intent(out) :: forecast(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: x_last, x_ar1
    call garch_midas_multi_step_forecast(param, spec, last_short_run_variance, &
      last_long_run_variance, horizon, forecast, status, x_last, x_ar1)
  end subroutine multi_step_ahead_pred

end module rumidas_forecast
