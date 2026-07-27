! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask_simulation
  use iso_fortran_env, only: int64
  use bidask_kinds, only: dp
  use bidask_types, only: ohlc_data
  use bidask_rng, only: rng_state
  implicit none
  private
  public :: simulate_ohlc, sim

  interface sim
    module procedure simulate_ohlc
  end interface sim

contains

  function simulate_ohlc(n, trades, probability, spread, volatility, overnight, &
      drift, signed_prices, seed) result(data)
    integer, intent(in), optional :: n, trades
    real(dp), intent(in), optional :: probability, spread, volatility, overnight, drift
    logical, intent(in), optional :: signed_prices
    integer(int64), intent(in), optional :: seed
    type(ohlc_data) :: data
    type(rng_state) :: rng
    real(dp) :: prob_value, spread_value, vol_value, overnight_value, drift_value
    real(dp) :: cumulative, z, price, previous, scale, mean_return
    real(dp), allocatable :: prices(:)
    logical, allocatable :: keep(:)
    logical :: use_sign
    integer :: n_value, trades_value, m, i, j, k, first, last, max_idx, min_idx
    real(dp) :: max_abs, min_abs

    n_value = 10000
    trades_value = 390
    if (present(n)) n_value = n
    if (present(trades)) trades_value = trades
    prob_value = 1.0_dp
    spread_value = 0.01_dp
    vol_value = 0.03_dp
    overnight_value = 0.0_dp
    drift_value = 0.0_dp
    use_sign = .false.
    if (present(probability)) prob_value = probability
    if (present(spread)) spread_value = spread
    if (present(volatility)) vol_value = volatility
    if (present(overnight)) overnight_value = overnight
    if (present(drift)) drift_value = drift
    if (present(signed_prices)) use_sign = signed_prices

    allocate(data%open(max(n_value, 0)), data%high(max(n_value, 0)), &
      data%low(max(n_value, 0)), data%close(max(n_value, 0)))
    if (n_value <= 0 .or. trades_value <= 0) return
    if (prob_value < 0.0_dp .or. prob_value > 1.0_dp) return
    if (spread_value < 0.0_dp .or. spread_value >= 2.0_dp) return
    if (vol_value < 0.0_dp .or. overnight_value < 0.0_dp) return

    if (present(seed)) then
      call rng%seed(seed)
    else
      call rng%seed(88172645463393265_int64)
    end if
    m = n_value * trades_value
    allocate(prices(m), keep(m))
    cumulative = 0.0_dp
    mean_return = drift_value / real(trades_value, dp)
    scale = vol_value / sqrt(real(trades_value, dp))
    do k = 1, m
      price = mean_return + scale * rng%normal()
      if (mod(k - 1, trades_value) == 0) price = price + overnight_value * rng%normal()
      cumulative = cumulative + price
      if (rng%bernoulli(0.5_dp)) then
        z = spread_value / 2.0_dp
      else
        z = -spread_value / 2.0_dp
      end if
      prices(k) = exp(cumulative) * (1.0_dp + z)
      if (use_sign) then
        if (z < 0.0_dp) then
          prices(k) = -prices(k)
        else if (spread_value <= 0.0_dp) then
          prices(k) = 0.0_dp
        end if
      end if
      keep(k) = rng%bernoulli(prob_value)
    end do

    previous = prices(1)
    do i = 1, n_value
      first = (i - 1) * trades_value + 1
      last = i * trades_value
      j = first
      do while (j <= last)
        if (keep(j)) exit
        j = j + 1
      end do
      if (j > last) then
        data%open(i) = previous
        data%high(i) = previous
        data%low(i) = previous
        data%close(i) = previous
        cycle
      end if
      data%open(i) = prices(j)
      max_idx = j
      min_idx = j
      max_abs = abs(prices(j))
      min_abs = abs(prices(j))
      data%close(i) = prices(j)
      do k = j + 1, last
        if (.not. keep(k)) cycle
        if (abs(prices(k)) > max_abs) then
          max_abs = abs(prices(k))
          max_idx = k
        end if
        if (abs(prices(k)) < min_abs) then
          min_abs = abs(prices(k))
          min_idx = k
        end if
        data%close(i) = prices(k)
      end do
      data%high(i) = prices(max_idx)
      data%low(i) = prices(min_idx)
      previous = data%close(i)
    end do
  end function simulate_ohlc

end module bidask_simulation
