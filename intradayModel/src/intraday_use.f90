! SPDX-License-Identifier: Apache-2.0
module intraday_use
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use intraday_kinds, only : dp
  use intraday_types, only : volume_model, volume_decomposition, kalman_output, &
                             intraday_ok, intraday_invalid_input, intraday_not_converged, &
                             intraday_numerical_failure
  use intraday_utils, only : clean_volume_data, parameters_valid, compute_error_metrics
  use intraday_kalman, only : uniss_kalman
  implicit none
  private

  public :: decompose_volume, forecast_volume

contains

  subroutine decompose_volume(purpose, model, data, result, burn_in_days)
    character(len=*), intent(in) :: purpose
    type(volume_model), intent(in) :: model
    real(dp), intent(in) :: data(:, :)
    type(volume_decomposition), intent(out) :: result
    integer, intent(in), optional :: burn_in_days

    real(dp), allocatable :: clean(:, :), log_data(:, :), original_all(:)
    real(dp), allocatable :: daily_all(:), dynamic_all(:), seasonal_all(:)
    type(kalman_output) :: kf
    integer :: status, burn, first, n, t, bin
    character(len=160) :: message
    character(len=:), allocatable :: mode
    logical :: is_forecast

    result%status = intraday_ok
    result%message = ''
    mode = lower_ascii(trim(adjustl(purpose)))
    is_forecast = mode == 'forecast'
    if (.not. is_forecast .and. mode /= 'analysis' .and. mode /= 'smooth') then
      result%status = intraday_invalid_input
      result%message = 'purpose must be analysis, smooth, or forecast'
      return
    end if
    result%is_forecast = is_forecast

    if (.not. model%converged) then
      result%status = intraday_not_converged
      result%message = 'all model parameters must be fitted or fixed before decomposition'
      return
    end if

    call clean_volume_data(data, clean, status=status, message=message)
    if (status /= intraday_ok) then
      result%status = status
      result%message = trim(message)
      return
    end if
    if (.not. parameters_valid(model%par, size(clean, 1))) then
      result%status = intraday_invalid_input
      result%message = 'model parameters do not match the data bins'
      return
    end if

    burn = 0
    if (present(burn_in_days)) burn = burn_in_days
    if (burn < 0 .or. burn > size(clean, 2)) then
      result%status = intraday_invalid_input
      result%message = 'burn_in_days must be between zero and the number of days'
      return
    end if
    if (.not. is_forecast) burn = 0

    allocate(log_data(size(clean, 1), size(clean, 2)))
    log_data = log(clean)
    call uniss_kalman(log_data, model%par, kf, smooth=.not. is_forecast)
    if (kf%status /= intraday_ok) then
      result%status = kf%status
      result%message = trim(kf%message)
      return
    end if

    n = size(clean)
    allocate(original_all(n), daily_all(n), dynamic_all(n), seasonal_all(n))
    original_all = reshape(clean, [n])
    do t = 1, n
      bin = mod(t - 1, size(clean, 1)) + 1
      seasonal_all(t) = exp(model%par%phi(bin))
      if (is_forecast) then
        daily_all(t) = exp(kf%x_pred(1, t))
        dynamic_all(t) = exp(kf%x_pred(2, t))
      else
        daily_all(t) = exp(kf%x_smooth(1, t))
        dynamic_all(t) = exp(kf%x_smooth(2, t))
      end if
    end do

    first = burn * size(clean, 1) + 1
    call copy_tail(original_all, first, result%original_signal)
    call copy_tail(daily_all, first, result%daily)
    call copy_tail(dynamic_all, first, result%dynamic)
    call copy_tail(seasonal_all, first, result%seasonal)
    allocate(result%fitted_signal(size(result%original_signal)))
    allocate(result%residual(size(result%original_signal)))
    result%fitted_signal = result%daily * result%dynamic * result%seasonal
    if (.not. all(ieee_is_finite(result%fitted_signal))) then
      result%status = intraday_numerical_failure
      result%message = 'component exponentiation overflowed'
      return
    end if
    result%residual = result%original_signal / result%fitted_signal
    result%error = compute_error_metrics(result%original_signal, result%fitted_signal)
  end subroutine decompose_volume

  subroutine forecast_volume(model, data, result, burn_in_days)
    type(volume_model), intent(in) :: model
    real(dp), intent(in) :: data(:, :)
    type(volume_decomposition), intent(out) :: result
    integer, intent(in), optional :: burn_in_days

    if (present(burn_in_days)) then
      call decompose_volume('forecast', model, data, result, burn_in_days)
    else
      call decompose_volume('forecast', model, data, result)
    end if
  end subroutine forecast_volume

  subroutine copy_tail(source, first, target)
    real(dp), intent(in) :: source(:)
    integer, intent(in) :: first
    real(dp), allocatable, intent(out) :: target(:)
    integer :: n

    n = max(0, size(source) - first + 1)
    allocate(target(n))
    if (n > 0) target = source(first:)
  end subroutine copy_tail

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code

    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module intraday_use
