! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_data
  use highfrequency_kinds, only: dp
  use highfrequency_types, only: liquidity_result
  implicit none
  private
  public :: make_returns, aggregate_last, aggregate_sum, previous_tick
  public :: refresh_time_pair, trade_direction, liquidity_measures
  public :: merge_same_timestamp, make_ohlcv

  interface make_returns
    module procedure make_returns_vector
    module procedure make_returns_matrix
  end interface make_returns

contains

  pure function make_returns_vector(prices, log_return) result(returns)
    real(dp), intent(in) :: prices(:)
    logical, intent(in), optional :: log_return
    real(dp), allocatable :: returns(:)
    logical :: use_log
    integer :: i, n
    n = size(prices)
    allocate(returns(max(0,n-1)))
    use_log = .true.
    if (present(log_return)) use_log = log_return
    do i = 1, n-1
      if (use_log) then
        if (prices(i) <= 0.0_dp .or. prices(i+1) <= 0.0_dp) then
          returns(i) = 0.0_dp
        else
          returns(i) = log(prices(i+1)/prices(i))
        end if
      else
        if (abs(prices(i)) <= tiny(1.0_dp)) then
          returns(i) = 0.0_dp
        else
          returns(i) = prices(i+1)/prices(i)-1.0_dp
        end if
      end if
    end do
  end function make_returns_vector

  pure function make_returns_matrix(prices, log_return) result(returns)
    real(dp), intent(in) :: prices(:,:)
    logical, intent(in), optional :: log_return
    real(dp), allocatable :: returns(:,:)
    logical :: use_log
    integer :: i, j, n, m
    n = size(prices,1)
    m = size(prices,2)
    allocate(returns(max(0,n-1),m))
    use_log = .true.
    if (present(log_return)) use_log = log_return
    do j = 1, m
      do i = 1, n-1
        if (use_log) then
          if (prices(i,j) <= 0.0_dp .or. prices(i+1,j) <= 0.0_dp) then
            returns(i,j) = 0.0_dp
          else
            returns(i,j) = log(prices(i+1,j)/prices(i,j))
          end if
        else
          if (abs(prices(i,j)) <= tiny(1.0_dp)) then
            returns(i,j) = 0.0_dp
          else
            returns(i,j) = prices(i+1,j)/prices(i,j)-1.0_dp
          end if
        end if
      end do
    end do
  end function make_returns_matrix

  subroutine aggregate_last(times, values, interval, out_times, out_values, nout)
    integer, intent(in) :: times(:), interval
    real(dp), intent(in) :: values(:,:)
    integer, allocatable, intent(out) :: out_times(:)
    real(dp), allocatable, intent(out) :: out_values(:,:)
    integer, intent(out) :: nout
    integer, allocatable :: tt(:)
    real(dp), allocatable :: vv(:,:)
    integer :: i, bin, previous_bin, n, m, count
    n = size(times)
    m = size(values,2)
    nout = 0
    if (n == 0 .or. interval <= 0 .or. size(values,1) /= n) then
      allocate(out_times(0), out_values(0,m))
      return
    end if
    allocate(tt(n), vv(n,m))
    previous_bin = floor(real(times(1),dp)/real(interval,dp))
    count = 1
    tt(count) = (previous_bin+1)*interval
    vv(count,:) = values(1,:)
    do i = 2, n
      bin = floor(real(times(i),dp)/real(interval,dp))
      if (bin == previous_bin) then
        vv(count,:) = values(i,:)
      else
        count = count + 1
        previous_bin = bin
        tt(count) = (bin+1)*interval
        vv(count,:) = values(i,:)
      end if
    end do
    nout = count
    allocate(out_times(count), out_values(count,m))
    out_times = tt(:count)
    out_values = vv(:count,:)
  end subroutine aggregate_last

  subroutine aggregate_sum(times, values, interval, out_times, out_values, nout)
    integer, intent(in) :: times(:), interval
    real(dp), intent(in) :: values(:,:)
    integer, allocatable, intent(out) :: out_times(:)
    real(dp), allocatable, intent(out) :: out_values(:,:)
    integer, intent(out) :: nout
    integer, allocatable :: tt(:)
    real(dp), allocatable :: vv(:,:)
    integer :: i, bin, previous_bin, n, m, count
    n = size(times)
    m = size(values,2)
    nout = 0
    if (n == 0 .or. interval <= 0 .or. size(values,1) /= n) then
      allocate(out_times(0), out_values(0,m))
      return
    end if
    allocate(tt(n), vv(n,m))
    previous_bin = floor(real(times(1),dp)/real(interval,dp))
    count = 1
    tt(count) = (previous_bin+1)*interval
    vv(count,:) = values(1,:)
    do i = 2, n
      bin = floor(real(times(i),dp)/real(interval,dp))
      if (bin == previous_bin) then
        vv(count,:) = vv(count,:) + values(i,:)
      else
        count = count + 1
        previous_bin = bin
        tt(count) = (bin+1)*interval
        vv(count,:) = values(i,:)
      end if
    end do
    nout = count
    allocate(out_times(count), out_values(count,m))
    out_times = tt(:count)
    out_values = vv(:count,:)
  end subroutine aggregate_sum

  subroutine previous_tick(source_times, source_values, target_times, target_values, available)
    integer, intent(in) :: source_times(:), target_times(:)
    real(dp), intent(in) :: source_values(:)
    real(dp), intent(out) :: target_values(size(target_times))
    logical, intent(out), optional :: available(size(target_times))
    integer :: i, j
    target_values = 0.0_dp
    if (present(available)) available = .false.
    j = 0
    do i = 1, size(target_times)
      do while (j < size(source_times))
        if (source_times(j+1) > target_times(i)) exit
        j = j + 1
      end do
      if (j > 0) then
        target_values(i) = source_values(j)
        if (present(available)) available(i) = .true.
      end if
    end do
  end subroutine previous_tick

  subroutine refresh_time_pair(t1, p1, t2, p2, times, x1, x2, nout)
    integer, intent(in) :: t1(:), t2(:)
    real(dp), intent(in) :: p1(:), p2(:)
    integer, allocatable, intent(out) :: times(:)
    real(dp), allocatable, intent(out) :: x1(:), x2(:)
    integer, intent(out) :: nout
    integer, allocatable :: tt(:)
    real(dp), allocatable :: xx(:), yy(:)
    integer :: i, j, k, rt, nmax
    nmax = min(size(t1), size(t2))
    allocate(tt(nmax), xx(nmax), yy(nmax))
    i = 1
    j = 1
    k = 0
    do while (i <= size(t1) .and. j <= size(t2))
      rt = max(t1(i), t2(j))
      do while (i < size(t1))
        if (t1(i+1) > rt) exit
        i = i + 1
      end do
      do while (j < size(t2))
        if (t2(j+1) > rt) exit
        j = j + 1
      end do
      if (t1(i) <= rt .and. t2(j) <= rt) then
        k = k + 1
        tt(k) = rt
        xx(k) = p1(i)
        yy(k) = p2(j)
      end if
      i = i + 1
      j = j + 1
    end do
    nout = k
    allocate(times(k), x1(k), x2(k))
    if (k > 0) then
      times = tt(:k)
      x1 = xx(:k)
      x2 = yy(:k)
    end if
  end subroutine refresh_time_pair

  subroutine trade_direction(price, bid, offer, direction)
    real(dp), intent(in) :: price(:), bid(:), offer(:)
    integer, intent(out) :: direction(size(price))
    real(dp) :: midpoint, change
    integer :: i, tick
    if (size(price) == 0) return
    tick = 1
    direction(1) = tick
    midpoint = 0.5_dp*(bid(1)+offer(1))
    if (price(1) < midpoint) direction(1) = -1
    if (price(1) > midpoint) direction(1) = 1
    tick = direction(1)
    do i = 2, size(price)
      change = price(i)-price(i-1)
      if (change > 0.0_dp) tick = 1
      if (change < 0.0_dp) tick = -1
      midpoint = 0.5_dp*(bid(i)+offer(i))
      if (price(i) < midpoint) then
        direction(i) = -1
      else if (price(i) > midpoint) then
        direction(i) = 1
      else
        direction(i) = tick
      end if
    end do
  end subroutine trade_direction

  subroutine liquidity_measures(price, size_trade, bid, offer, bid_size, offer_size, window, result)
    real(dp), intent(in) :: price(:), size_trade(:), bid(:), offer(:), bid_size(:), offer_size(:)
    integer, intent(in) :: window
    type(liquidity_result), intent(out) :: result
    integer :: n, i, j
    real(dp) :: denom
    n = size(price)
    allocate(result%direction(n), result%midpoint(n), result%effective_spread(n))
    allocate(result%realized_spread(n), result%price_impact(n))
    allocate(result%proportional_effective_spread(n), result%proportional_realized_spread(n))
    allocate(result%proportional_price_impact(n), result%depth_imbalance_difference(n))
    allocate(result%depth_imbalance_ratio(n), result%signed_trade_size(n))
    allocate(result%quoted_spread(n), result%proportional_quoted_spread(n))
    call trade_direction(price,bid,offer,result%direction)
    result%midpoint = 0.5_dp*(bid+offer)
    do i = 1, n
      result%effective_spread(i) = 2.0_dp*real(result%direction(i),dp)*(price(i)-result%midpoint(i))
      j = i + max(0,window)
      if (j <= n) then
        result%realized_spread(i) = 2.0_dp*real(result%direction(i),dp)*(price(i)-result%midpoint(j))
      else
        result%realized_spread(i) = 0.0_dp
      end if
      result%price_impact(i) = 0.5_dp*(result%effective_spread(i)-result%realized_spread(i))
      if (abs(result%midpoint(i)) > tiny(1.0_dp)) then
        result%proportional_effective_spread(i) = result%effective_spread(i)/result%midpoint(i)
        result%proportional_realized_spread(i) = result%realized_spread(i)/result%midpoint(i)
        result%proportional_price_impact(i) = result%price_impact(i)/result%midpoint(i)
      else
        result%proportional_effective_spread(i) = 0.0_dp
        result%proportional_realized_spread(i) = 0.0_dp
        result%proportional_price_impact(i) = 0.0_dp
      end if
      denom = offer_size(i)+bid_size(i)
      if (abs(denom) > tiny(1.0_dp)) then
        result%depth_imbalance_difference(i) = real(result%direction(i),dp) * &
          (offer_size(i)-bid_size(i))/denom
      else
        result%depth_imbalance_difference(i) = 0.0_dp
      end if
      if (bid_size(i) > 0.0_dp .and. offer_size(i) > 0.0_dp) then
        if (result%direction(i) == 1) then
          result%depth_imbalance_ratio(i) = offer_size(i)/bid_size(i)
        else
          result%depth_imbalance_ratio(i) = bid_size(i)/offer_size(i)
        end if
      else
        result%depth_imbalance_ratio(i) = 0.0_dp
      end if
      result%signed_trade_size(i) = real(result%direction(i),dp)*size_trade(i)
      result%quoted_spread(i) = offer(i)-bid(i)
      if (abs(result%midpoint(i)) > tiny(1.0_dp)) then
        result%proportional_quoted_spread(i) = result%quoted_spread(i)/result%midpoint(i)
      else
        result%proportional_quoted_spread(i) = 0.0_dp
      end if
    end do
  end subroutine liquidity_measures

  subroutine merge_same_timestamp(times, values, out_times, out_values, nout, sum_values)
    integer, intent(in) :: times(:)
    real(dp), intent(in) :: values(:,:)
    integer, allocatable, intent(out) :: out_times(:)
    real(dp), allocatable, intent(out) :: out_values(:,:)
    integer, intent(out) :: nout
    logical, intent(in), optional :: sum_values
    integer, allocatable :: tt(:)
    real(dp), allocatable :: vv(:,:)
    logical :: use_sum
    integer :: i, count, n, m
    n = size(times)
    m = size(values,2)
    allocate(tt(n), vv(n,m))
    use_sum = .false.
    if (present(sum_values)) use_sum = sum_values
    count = 0
    do i = 1, n
      if (count == 0 .or. times(i) /= tt(count)) then
        count = count+1
        tt(count) = times(i)
        vv(count,:) = values(i,:)
      else if (use_sum) then
        vv(count,:) = vv(count,:)+values(i,:)
      else
        vv(count,:) = 0.5_dp*(vv(count,:)+values(i,:))
      end if
    end do
    nout = count
    allocate(out_times(count),out_values(count,m))
    if (count > 0) then
      out_times=tt(:count)
      out_values=vv(:count,:)
    end if
  end subroutine merge_same_timestamp

  subroutine make_ohlcv(times, price, volume, interval, out_times, ohlcv, nout)
    integer, intent(in) :: times(:), interval
    real(dp), intent(in) :: price(:), volume(:)
    integer, allocatable, intent(out) :: out_times(:)
    real(dp), allocatable, intent(out) :: ohlcv(:,:)
    integer, intent(out) :: nout
    integer, allocatable :: tt(:)
    real(dp), allocatable :: data(:,:)
    integer :: i, count, bin, previous_bin, n
    n=size(times)
    allocate(tt(n),data(n,5))
    count=0
    previous_bin=-huge(1)
    do i=1,n
      bin=floor(real(times(i),dp)/real(interval,dp))
      if (bin /= previous_bin) then
        count=count+1
        previous_bin=bin
        tt(count)=(bin+1)*interval
        data(count,1)=price(i)
        data(count,2)=price(i)
        data(count,3)=price(i)
        data(count,4)=price(i)
        data(count,5)=volume(i)
      else
        data(count,2)=max(data(count,2),price(i))
        data(count,3)=min(data(count,3),price(i))
        data(count,4)=price(i)
        data(count,5)=data(count,5)+volume(i)
      end if
    end do
    nout=count
    allocate(out_times(count),ohlcv(count,5))
    if(count>0)then
      out_times=tt(:count)
      ohlcv=data(:count,:)
    end if
  end subroutine make_ohlcv

end module highfrequency_data
