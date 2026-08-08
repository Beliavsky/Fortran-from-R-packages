module nfcp_forecast
  use nfcp_types, only : dp, nfcp_model_t, nfcp_ok, nfcp_invalid_input
  use nfcp_math, only : nfcp_a_t, nfcp_covariance, nfcp_seasonality, normal_quantile
  implicit none
  private
  public :: spot_price_forecast, futures_price_forecast

contains

  subroutine spot_price_forecast(model, initial_state, times, forecast, percentiles, seasonal_trend, status)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: initial_state(:), times(:)
    real(dp), allocatable, intent(out) :: forecast(:,:)
    real(dp), intent(in), optional :: percentiles(:)
    real(dp), intent(in), optional :: seasonal_trend
    integer, intent(out), optional :: status
    integer :: i, j, nq
    real(dp) :: trend, mean_log, variance
    real(dp), allocatable :: q(:), covariance(:,:)

    if (present(status)) status = nfcp_invalid_input
    if (size(initial_state) /= model%n_factors .or. any(times < 0.0_dp)) then
      allocate(forecast(0,0))
      return
    end if
    trend = 0.0_dp
    if (present(seasonal_trend)) trend = seasonal_trend
    nq = 1
    if (present(percentiles)) nq = 1 + size(percentiles)
    allocate(q(nq), forecast(size(times),nq), covariance(model%n_factors,model%n_factors))
    q(1) = 0.5_dp
    if (present(percentiles)) then
      if (any(percentiles < 0.0_dp) .or. any(percentiles > 1.0_dp)) then
        deallocate(forecast); allocate(forecast(0,0)); return
      end if
      q(2:) = percentiles
    end if

    do i = 1, size(times)
      mean_log = model%equilibrium + model%mu*times(i) + nfcp_seasonality(model,times(i)+trend)
      do j = 1, model%n_factors
        mean_log = mean_log + initial_state(j)*exp(-model%kappa(j)*times(i))
      end do
      call nfcp_covariance(model,times(i),covariance)
      variance = sum(covariance)
      mean_log = mean_log + 0.5_dp*variance
      do j = 1, nq
        forecast(i,j) = exp(mean_log + sqrt(max(0.0_dp,variance))*normal_quantile(q(j)))
      end do
    end do
    if (present(status)) status = nfcp_ok
  end subroutine spot_price_forecast

  subroutine futures_price_forecast(model, initial_state, observation_time, maturities, forecast, &
                                    percentiles, seasonal_trend, status)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: initial_state(:), observation_time, maturities(:)
    real(dp), allocatable, intent(out) :: forecast(:,:)
    real(dp), intent(in), optional :: percentiles(:)
    real(dp), intent(in), optional :: seasonal_trend
    integer, intent(out), optional :: status
    integer :: i, j, k, nq
    real(dp) :: trend, mean_log, variance, ks
    real(dp), allocatable :: q(:), covariance(:,:)

    if (present(status)) status = nfcp_invalid_input
    if (size(initial_state) /= model%n_factors .or. observation_time < 0.0_dp .or. &
        any(maturities < observation_time)) then
      allocate(forecast(0,0)); return
    end if
    trend = 0.0_dp
    if (present(seasonal_trend)) trend = seasonal_trend
    nq = 1
    if (present(percentiles)) nq = 1+size(percentiles)
    allocate(q(nq), forecast(size(maturities),nq), covariance(model%n_factors,model%n_factors))
    q(1)=0.5_dp
    if (present(percentiles)) then
      if (any(percentiles < 0.0_dp) .or. any(percentiles > 1.0_dp)) then
        deallocate(forecast); allocate(forecast(0,0)); return
      end if
      q(2:)=percentiles
    end if

    call nfcp_covariance(model,observation_time,covariance)
    do i = 1, size(maturities)
      mean_log = model%equilibrium + model%mu_rn*observation_time + &
                 nfcp_seasonality(model,maturities(i)+trend) + &
                 nfcp_a_t(model,maturities(i)-observation_time)
      do j = 1, model%n_factors
        mean_log = mean_log + initial_state(j)*exp(-model%kappa(j)*maturities(i))
      end do
      variance = 0.0_dp
      do j = 1, model%n_factors
        do k = 1, model%n_factors
          ks = model%kappa(j)+model%kappa(k)
          variance = variance + covariance(j,k)*exp(-ks*(maturities(i)-observation_time))
        end do
      end do
      do j = 1, nq
        forecast(i,j)=exp(mean_log+sqrt(max(0.0_dp,variance))*normal_quantile(q(j)))
      end do
    end do
    if (present(status)) status=nfcp_ok
  end subroutine futures_price_forecast

end module nfcp_forecast
