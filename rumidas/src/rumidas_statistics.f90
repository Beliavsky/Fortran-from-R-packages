! SPDX-License-Identifier: GPL-3.0-only
module rumidas_statistics
  use rumidas_kinds, only: dp
  implicit none
  private
  public :: information_criteria, volatility_loss_functions, sample_mean, sample_variance

contains

  pure function sample_mean(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function sample_mean

  pure function sample_variance(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, mean_x
    if (size(x) <= 1) then
      value = 0.0_dp
    else
      mean_x = sample_mean(x)
      value = sum((x - mean_x) ** 2) / real(size(x) - 1, dp)
    end if
  end function sample_variance

  pure subroutine information_criteria(loglik, number_parameters, number_observations, aic, bic)
    real(dp), intent(in) :: loglik
    integer, intent(in) :: number_parameters, number_observations
    real(dp), intent(out) :: aic, bic
    aic = 2.0_dp * real(number_parameters, dp) - 2.0_dp * loglik
    bic = real(number_parameters, dp) * log(real(max(number_observations, 1), dp)) - 2.0_dp * loglik
  end subroutine information_criteria

  pure subroutine volatility_loss_functions(volatility_estimate, volatility_proxy, mse_percent, qlike)
    real(dp), intent(in) :: volatility_estimate(:), volatility_proxy(:)
    real(dp), intent(out) :: mse_percent, qlike
    integer :: n
    n = min(size(volatility_estimate), size(volatility_proxy))
    if (n <= 0 .or. any(volatility_estimate(1:n) <= 0.0_dp)) then
      mse_percent = huge(1.0_dp)
      qlike = huge(1.0_dp)
      return
    end if
    mse_percent = 100.0_dp * sum((volatility_estimate(1:n) - volatility_proxy(1:n)) ** 2) / real(n, dp)
    qlike = sum(log(volatility_estimate(1:n)) + volatility_proxy(1:n) / volatility_estimate(1:n)) / real(n, dp)
  end subroutine volatility_loss_functions

end module rumidas_statistics
