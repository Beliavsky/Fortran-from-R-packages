! SPDX-License-Identifier: MIT
module mfgarch_prediction
  use mfgarch_kinds, only : dp
  use mfgarch_components, only : forecast_garch
  use mfgarch_status, only : mfgarch_success, mfgarch_invalid_argument
  use mfgarch_types, only : mfgarch_model
  implicit none
  private

  public :: predict_variance

contains

  subroutine predict_variance(model, horizons, tau_forecast, last_return, conditional_g, &
      conditional_tau, forecasts, status)
    type(mfgarch_model), intent(in) :: model
    integer, intent(in) :: horizons(:)
    real(dp), intent(in) :: tau_forecast, last_return, conditional_g, conditional_tau
    real(dp), allocatable, intent(out) :: forecasts(:)
    integer, intent(out) :: status
    real(dp) :: standardized_return, omega, gamma_value
    integer :: i

    status = mfgarch_success
    if (any(horizons < 1) .or. tau_forecast <= 0.0_dp .or. conditional_g <= 0.0_dp .or. &
        conditional_tau <= 0.0_dp) then
      status = mfgarch_invalid_argument
      allocate(forecasts(0))
      return
    end if
    allocate(forecasts(size(horizons)))
    gamma_value = merge(model%gamma, 0.0_dp, model%asymmetric)
    omega = 1.0_dp - model%alpha - model%beta - 0.5_dp * gamma_value
    standardized_return = (last_return - model%mu) / sqrt(conditional_tau)
    do i = 1, size(horizons)
      forecasts(i) = tau_forecast * forecast_garch(omega, model%alpha, model%beta, &
        gamma_value, conditional_g, standardized_return, horizons(i))
    end do
  end subroutine predict_variance

end module mfgarch_prediction
