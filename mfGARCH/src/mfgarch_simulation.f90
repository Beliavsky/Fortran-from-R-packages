! SPDX-License-Identifier: MIT
module mfgarch_simulation
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch_kinds, only : dp
  use mfgarch_components, only : beta_weights
  use mfgarch_math, only : rolling_mean, finite_value
  use mfgarch_random, only : mfgarch_rng
  use mfgarch_status, only : mfgarch_success, mfgarch_invalid_argument, &
    mfgarch_numerical_error
  use mfgarch_types, only : mfgarch_model, mfgarch_simulation_type => mfgarch_simulation
  implicit none
  private

  public :: simulate_mfgarch, simulate_mfgarch_rv_dependent
  public :: simulate_mfgarch_diffusion

contains

  subroutine simulate_mfgarch(n_days, model, psi, sigma_psi, low_frequency, n_intraday, &
      simulation, status, seed, student_t_df, correlation)
    integer, intent(in) :: n_days, low_frequency, n_intraday
    type(mfgarch_model), intent(in) :: model
    real(dp), intent(in) :: psi, sigma_psi
    type(mfgarch_simulation_type), intent(out) :: simulation
    integer, intent(out) :: status
    integer(int64), intent(in), optional :: seed
    real(dp), intent(in), optional :: student_t_df, correlation
    type(mfgarch_rng) :: rng
    real(dp), allocatable :: innovations(:), standardized_intraday(:), daily_standardized(:)
    real(dp), allocatable :: g_full(:), x_low(:), tau_low(:), tau_full(:), ret_full(:)
    real(dp), allocatable :: block_shock(:), weights(:)
    real(dp) :: corr, z, scale_t, omega, shock, log_tau
    integer :: burnin, total_days, total_intraday, nlow, day, j, p, first, last

    status = mfgarch_success
    if (n_days <= 0 .or. low_frequency <= 0 .or. n_intraday <= 0 .or. &
        mod(n_intraday,48) /= 0 .or. mod(n_days,low_frequency) /= 0 .or. &
        abs(psi) >= 1.0_dp .or. sigma_psi < 0.0_dp .or. model%k < 0) then
      status = mfgarch_invalid_argument
      return
    end if
    if (present(student_t_df)) then
      if (student_t_df <= 2.0_dp) then
        status = mfgarch_invalid_argument
        return
      end if
    end if
    corr = 0.0_dp
    if (present(correlation)) corr = correlation
    if (abs(corr) > 1.0_dp) then
      status = mfgarch_invalid_argument
      return
    end if
    if (present(seed)) then
      call rng%seed(seed)
    else
      call rng%seed(2718281828459_int64)
    end if

    burnin = 3 * low_frequency * model%k
    total_days = n_days + burnin
    total_intraday = total_days * n_intraday
    nlow = total_days / low_frequency
    allocate(innovations(total_intraday), standardized_intraday(total_intraday), &
      daily_standardized(total_days), g_full(total_days), x_low(nlow), tau_low(nlow), &
      tau_full(total_days), ret_full(total_intraday), block_shock(nlow))

    scale_t = 1.0_dp
    if (present(student_t_df)) scale_t = sqrt(student_t_df / (student_t_df - 2.0_dp))
    do j = 1, total_intraday
      if (present(student_t_df)) then
        innovations(j) = rng%student_t(student_t_df) / scale_t
      else
        innovations(j) = rng%normal()
      end if
    end do

    omega = 1.0_dp - model%alpha - model%beta - 0.5_dp * merge(model%gamma, 0.0_dp, model%asymmetric)
    if (omega <= 0.0_dp) then
      status = mfgarch_invalid_argument
      return
    end if
    g_full(1) = 1.0_dp
    do day = 1, total_days
      first = (day - 1) * n_intraday + 1
      last = day * n_intraday
      standardized_intraday(first:last) = innovations(first:last) * sqrt(g_full(day) / real(n_intraday, dp))
      daily_standardized(day) = sum(standardized_intraday(first:last))
      if (day < total_days) then
        shock = model%alpha * daily_standardized(day)**2
        if (model%asymmetric .and. daily_standardized(day) < 0.0_dp) &
          shock = shock + model%gamma * daily_standardized(day)**2
        g_full(day+1) = omega + shock + model%beta * g_full(day)
        if (.not. finite_value(g_full(day+1)) .or. g_full(day+1) <= 0.0_dp) then
          status = mfgarch_numerical_error
          return
        end if
      end if
    end do

    do p = 1, nlow
      first = (p - 1) * low_frequency + 1
      last = p * low_frequency
      block_shock(p) = sum(daily_standardized(first:last)) / sqrt(real(low_frequency, dp))
    end do
    x_low = 0.0_dp
    do p = 2, nlow
      z = corr * block_shock(p) + sqrt(max(0.0_dp, 1.0_dp - corr*corr)) * rng%normal()
      x_low(p) = psi * x_low(p-1) + sigma_psi * z
    end do

    if (model%k == 0) then
      tau_low = exp(model%m)
    else
      call beta_weights(model%k, model%w1, model%w2, weights, status)
      if (status /= mfgarch_success) return
      tau_low = exp(model%m)
      do p = model%k + 1, nlow
        log_tau = model%m
        do j = 1, model%k
          log_tau = log_tau + model%theta * weights(j) * x_low(p-j)
        end do
        tau_low(p) = exp(log_tau)
      end do
    end if
    do day = 1, total_days
      tau_full(day) = tau_low((day - 1) / low_frequency + 1)
      first = (day - 1) * n_intraday + 1
      last = day * n_intraday
      ret_full(first:last) = standardized_intraday(first:last) * sqrt(tau_full(day)) + &
        model%mu / real(n_intraday, dp)
    end do

    call finalize_simulation(ret_full, tau_full, g_full, x_low, burnin, n_days, &
      low_frequency, n_intraday, simulation, status)
  end subroutine simulate_mfgarch

  subroutine simulate_mfgarch_rv_dependent(n_days, model, low_frequency, n_intraday, &
      use_realized_volatility, simulation, status, seed)
    integer, intent(in) :: n_days, low_frequency, n_intraday
    type(mfgarch_model), intent(in) :: model
    logical, intent(in) :: use_realized_volatility
    type(mfgarch_simulation_type), intent(out) :: simulation
    integer, intent(out) :: status
    integer(int64), intent(in), optional :: seed
    type(mfgarch_rng) :: rng
    real(dp), allocatable :: ret_full(:), tau_full(:), g_full(:), x_low(:), tau_low(:), weights(:)
    real(dp) :: omega, shock, log_tau, standardized_previous
    integer :: burnin, total_days, total_intraday, nlow, day, p, j, first, last

    status = mfgarch_success
    if (n_days <= 0 .or. low_frequency <= 0 .or. n_intraday <= 0 .or. &
        mod(n_intraday,48) /= 0 .or. mod(n_days,low_frequency) /= 0 .or. model%k < 1) then
      status = mfgarch_invalid_argument
      return
    end if
    if (present(seed)) then
      call rng%seed(seed)
    else
      call rng%seed(3141592653589_int64)
    end if
    burnin = 10 * low_frequency * model%k
    total_days = n_days + burnin
    total_intraday = total_days * n_intraday
    nlow = total_days / low_frequency
    allocate(ret_full(total_intraday), tau_full(total_days), g_full(total_days), &
      x_low(nlow), tau_low(nlow))
    call beta_weights(model%k, model%w1, model%w2, weights, status)
    if (status /= mfgarch_success) return

    omega = 1.0_dp - model%alpha - model%beta - 0.5_dp * merge(model%gamma, 0.0_dp, model%asymmetric)
    if (omega <= 0.0_dp) then
      status = mfgarch_invalid_argument
      return
    end if
    x_low = 1.0_dp
    tau_low = exp(model%m)
    g_full(1) = 1.0_dp

    do day = 1, total_days
      p = (day - 1) / low_frequency + 1
      tau_full(day) = tau_low(p)
      first = (day - 1) * n_intraday + 1
      last = day * n_intraday
      do j = first, last
        ret_full(j) = rng%normal() * sqrt(g_full(day) * tau_full(day) / real(n_intraday, dp)) + &
          model%mu / real(n_intraday, dp)
      end do
      if (day < total_days) then
        standardized_previous = (sum(ret_full(first:last)) - model%mu) / sqrt(tau_full(day))
        shock = model%alpha * standardized_previous**2
        if (model%asymmetric .and. standardized_previous < 0.0_dp) &
          shock = shock + model%gamma * standardized_previous**2
        g_full(day+1) = omega + shock + model%beta * g_full(day)
        if (g_full(day+1) <= 0.0_dp .or. .not. finite_value(g_full(day+1))) then
          status = mfgarch_numerical_error
          return
        end if
      end if
      if (mod(day, low_frequency) == 0) then
        p = day / low_frequency
        first = (day - low_frequency) * n_intraday + 1
        last = day * n_intraday
        x_low(p) = 0.0_dp
        do j = day - low_frequency + 1, day
          x_low(p) = x_low(p) + sum(ret_full((j-1)*n_intraday+1:j*n_intraday))**2
        end do
        if (use_realized_volatility) x_low(p) = sqrt(x_low(p) / real(low_frequency, dp))
        if (p < nlow) then
          log_tau = model%m
          do j = 1, min(model%k, p)
            log_tau = log_tau + model%theta * weights(j) * x_low(p+1-j)
          end do
          tau_low(p+1) = exp(log_tau)
        end if
      end if
    end do

    call finalize_simulation(ret_full, tau_full, g_full, x_low, burnin, n_days, &
      low_frequency, n_intraday, simulation, status)
  end subroutine simulate_mfgarch_rv_dependent

  subroutine simulate_mfgarch_diffusion(n_days, model, psi, sigma_psi, low_frequency, &
      n_intraday, simulation, status, seed)
    integer, intent(in) :: n_days, low_frequency, n_intraday
    type(mfgarch_model), intent(in) :: model
    real(dp), intent(in) :: psi, sigma_psi
    type(mfgarch_simulation_type), intent(out) :: simulation
    integer, intent(out) :: status
    integer(int64), intent(in), optional :: seed
    type(mfgarch_rng) :: rng
    real(dp), allocatable :: raw_intraday(:), ret_full(:), tau_full(:), g_full(:)
    real(dp), allocatable :: x_low(:), tau_low(:), weights(:)
    real(dp) :: theta_short, lambda, h, price, previous_sample, z, log_tau, persistence
    integer :: burnin, total_days, nlow, sampling, delta, total_trades
    integer :: trade, day, p, j, intraday_index, first, last

    status = mfgarch_success
    persistence = model%alpha + model%beta
    if (n_days <= 0 .or. low_frequency <= 0 .or. n_intraday <= 0 .or. &
        mod(n_intraday,6) /= 0 .or. mod(n_days,low_frequency) /= 0 .or. &
        model%k < 0 .or. persistence <= 0.0_dp .or. persistence >= 1.0_dp .or. &
        model%alpha <= 0.0_dp .or. abs(psi) >= 1.0_dp .or. sigma_psi < 0.0_dp) then
      status = mfgarch_invalid_argument
      return
    end if
    if (present(seed)) then
      call rng%seed(seed)
    else
      call rng%seed(1618033988749_int64)
    end if
    burnin = 2 * low_frequency * model%k
    total_days = n_days + burnin
    nlow = total_days / low_frequency
    sampling = 20
    delta = n_intraday * sampling
    total_trades = total_days * delta
    allocate(raw_intraday(total_days*n_intraday), ret_full(total_days*n_intraday), &
      tau_full(total_days), g_full(total_days), x_low(nlow), tau_low(nlow))

    theta_short = -log(persistence)
    lambda = 2.0_dp * log(persistence)**2 / &
      ((((1.0_dp - persistence**2) * (1.0_dp - model%beta)**2) / &
      (model%alpha * (1.0_dp - model%beta * persistence))) + &
      6.0_dp * log(persistence) + 2.0_dp * log(persistence)**2 + &
      4.0_dp * (1.0_dp - persistence))
    if (.not. finite_value(lambda) .or. lambda <= 0.0_dp) then
      status = mfgarch_numerical_error
      return
    end if

    h = 0.1_dp
    price = 0.0_dp
    previous_sample = 0.0_dp
    intraday_index = 0
    g_full = 0.0_dp
    do trade = 1, total_trades
      if (trade > 1) then
        z = rng%normal()
        h = theta_short / real(delta, dp) + h * (1.0_dp - theta_short / real(delta, dp) + &
          sqrt(2.0_dp * lambda * theta_short / real(delta, dp)) * z)
        h = max(h, 1.0e-12_dp)
      end if
      price = price + sqrt(h / real(delta, dp)) * rng%normal()
      day = (trade - 1) / delta + 1
      g_full(day) = g_full(day) + h / real(delta, dp)
      if (mod(trade, sampling) == 0) then
        intraday_index = intraday_index + 1
        raw_intraday(intraday_index) = price - previous_sample
        previous_sample = price
      end if
    end do

    x_low = 0.0_dp
    do p = 2, nlow
      x_low(p) = psi * x_low(p-1) + sigma_psi * rng%normal()
    end do
    if (model%k == 0) then
      tau_low = exp(model%m)
    else
      call beta_weights(model%k, model%w1, model%w2, weights, status)
      if (status /= mfgarch_success) return
      tau_low = exp(model%m)
      do p = model%k + 1, nlow
        log_tau = model%m
        do j = 1, model%k
          log_tau = log_tau + model%theta * weights(j) * x_low(p-j)
        end do
        tau_low(p) = exp(log_tau)
      end do
    end if
    do day = 1, total_days
      p = (day - 1) / low_frequency + 1
      tau_full(day) = tau_low(p)
      first = (day - 1) * n_intraday + 1
      last = day * n_intraday
      ret_full(first:last) = raw_intraday(first:last) * sqrt(tau_full(day)) + &
        model%mu / real(n_intraday, dp)
    end do

    call finalize_simulation(ret_full, tau_full, g_full, x_low, burnin, n_days, &
      low_frequency, n_intraday, simulation, status, half_hour_groups=n_intraday/6)
  end subroutine simulate_mfgarch_diffusion

  subroutine finalize_simulation(ret_full, tau_full, g_full, x_low, burnin, n_days, &
      low_frequency, n_intraday, simulation, status, half_hour_groups)
    real(dp), intent(in) :: ret_full(:), tau_full(:), g_full(:), x_low(:)
    integer, intent(in) :: burnin, n_days, low_frequency, n_intraday
    type(mfgarch_simulation_type), intent(out) :: simulation
    integer, intent(out) :: status
    integer, intent(in), optional :: half_hour_groups
    integer :: day, source_day, first, last, group, group_size, ngroups, p, offset
    real(dp) :: group_return
    real(dp), allocatable :: full_daily_returns(:), full_rv(:), full_rv_half(:)

    status = mfgarch_success
    allocate(full_daily_returns(size(tau_full)), full_rv(size(tau_full)), full_rv_half(size(tau_full)))
    ngroups = 48
    if (present(half_hour_groups)) ngroups = half_hour_groups
    if (ngroups <= 0 .or. mod(n_intraday,ngroups) /= 0) then
      status = mfgarch_invalid_argument
      return
    end if
    group_size = n_intraday / ngroups
    do day = 1, size(tau_full)
      first = (day - 1) * n_intraday + 1
      last = day * n_intraday
      full_daily_returns(day) = sum(ret_full(first:last))
      full_rv(day) = sum(ret_full(first:last)**2)
      full_rv_half(day) = 0.0_dp
      do group = 1, ngroups
        offset = first + (group - 1) * group_size
        group_return = sum(ret_full(offset:offset+group_size-1))
        full_rv_half(day) = full_rv_half(day) + group_return**2
      end do
    end do

    allocate(simulation%returns(n_days), simulation%tau(n_days), simulation%g(n_days), &
      simulation%covariate(n_days), simulation%realized_variance(n_days), &
      simulation%realized_variance_half_hour(n_days), simulation%low_frequency_period(n_days), &
      simulation%intraday_returns(n_days*n_intraday))
    do day = 1, n_days
      source_day = burnin + day
      simulation%returns(day) = full_daily_returns(source_day)
      simulation%tau(day) = tau_full(source_day)
      simulation%g(day) = g_full(source_day)
      p = (source_day - 1) / low_frequency + 1
      simulation%covariate(day) = x_low(p)
      simulation%low_frequency_period(day) = (day - 1) / low_frequency + 1
      simulation%realized_variance(day) = full_rv(source_day)
      simulation%realized_variance_half_hour(day) = full_rv_half(source_day)
      first = (source_day - 1) * n_intraday + 1
      last = source_day * n_intraday
      offset = (day - 1) * n_intraday + 1
      simulation%intraday_returns(offset:offset+n_intraday-1) = ret_full(first:last)
    end do
    call rolling_mean(simulation%realized_variance, 5, simulation%rv_5, status)
    if (status /= mfgarch_success) return
    call rolling_mean(simulation%realized_variance, 22, simulation%rv_22, status)
    if (status /= mfgarch_success) return
    call rolling_mean(simulation%realized_variance_half_hour, 5, simulation%rv_half_hour_5, status)
    if (status /= mfgarch_success) return
    call rolling_mean(simulation%realized_variance_half_hour, 22, simulation%rv_half_hour_22, status)
  end subroutine finalize_simulation

end module mfgarch_simulation
