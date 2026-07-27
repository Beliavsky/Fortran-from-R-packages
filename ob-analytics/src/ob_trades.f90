! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_trades
   use ob_kinds, only : dp, i8
   use ob_types, only : event_t, trade_t, impact_t, side_bid, side_ask, trade_buy, trade_sell
   use ob_utils, only : weighted_average
   implicit none
   private
   public :: match_trades, trade_impacts

contains

   function match_trades(events, correct_large_jumps, jump_threshold) result(trades)
      type(event_t), intent(in) :: events(:)
      logical, intent(in), optional :: correct_large_jumps
      real(dp), intent(in), optional :: jump_threshold
      type(trade_t), allocatable :: trades(:)
      logical :: correct
      real(dp) :: threshold
      integer :: i, j, n, ask_idx
      logical :: bid_maker

      correct = .true.
      threshold = 10.0_dp
      if (present(correct_large_jumps)) correct = correct_large_jumps
      if (present(jump_threshold)) threshold = jump_threshold
      n = count(events%side == side_bid .and. events%matching_event > 0)
      allocate(trades(n))
      n = 0
      do i = 1, size(events)
         if (events(i)%side /= side_bid .or. events(i)%matching_event <= 0) cycle
         ask_idx = find_event(events, events(i)%matching_event)
         if (ask_idx == 0) cycle
         if (events(ask_idx)%side /= side_ask) cycle
         n = n + 1
         bid_maker = events(i)%exchange_timestamp_ms < events(ask_idx)%exchange_timestamp_ms .or. &
            (events(i)%exchange_timestamp_ms == events(ask_idx)%exchange_timestamp_ms .and. &
             events(i)%id < events(ask_idx)%id)
         trades(n)%timestamp_ms = min(events(i)%timestamp_ms, events(ask_idx)%timestamp_ms)
         trades(n)%price = merge(events(i)%price, events(ask_idx)%price, bid_maker)
         trades(n)%volume = events(i)%fill
         trades(n)%direction = merge(trade_sell, trade_buy, bid_maker)
         trades(n)%maker_event_id = merge(events(i)%event_id, events(ask_idx)%event_id, bid_maker)
         trades(n)%taker_event_id = merge(events(ask_idx)%event_id, events(i)%event_id, bid_maker)
         trades(n)%maker_id = merge(events(i)%id, events(ask_idx)%id, bid_maker)
         trades(n)%taker_id = merge(events(ask_idx)%id, events(i)%id, bid_maker)
      end do
      if (n < size(trades)) trades = trades(1:n)
      call sort_trades_time(trades)

      if (correct) then
         do i = 2, size(trades)
            if (abs(trades(i)%price-trades(i-1)%price) > threshold) then
               j = find_event(events, trades(i)%taker_event_id)
               if (j > 0) call swap_trade_roles(trades(i), events(j))
            end if
         end do
      end if
   end function match_trades

   function trade_impacts(trades) result(impacts)
      type(trade_t), intent(in) :: trades(:)
      type(impact_t), allocatable :: impacts(:)
      integer(i8), allocatable :: takers(:)
      real(dp), allocatable :: price(:), volume(:)
      integer :: i, j, n, k

      allocate(takers(size(trades)))
      n = 0
      do i = 1, size(trades)
         if (.not. any(takers(1:n) == trades(i)%taker_id)) then
            n = n + 1
            takers(n) = trades(i)%taker_id
         end if
      end do
      allocate(impacts(n))
      do i = 1, n
         k = count(trades%taker_id == takers(i))
         allocate(price(k), volume(k))
         k = 0
         do j = 1, size(trades)
            if (trades(j)%taker_id /= takers(i)) cycle
            k = k + 1
            price(k) = trades(j)%price
            volume(k) = trades(j)%volume
            if (k == 1) then
               impacts(i)%start_time_ms = trades(j)%timestamp_ms
               impacts(i)%end_time_ms = trades(j)%timestamp_ms
               impacts(i)%direction = trades(j)%direction
            else
               impacts(i)%start_time_ms = min(impacts(i)%start_time_ms, trades(j)%timestamp_ms)
               impacts(i)%end_time_ms = max(impacts(i)%end_time_ms, trades(j)%timestamp_ms)
            end if
         end do
         impacts(i)%id = takers(i)
         impacts(i)%min_price = minval(price)
         impacts(i)%max_price = maxval(price)
         impacts(i)%vwap = anint(100.0_dp*weighted_average(price,volume))/100.0_dp
         impacts(i)%hits = size(price)
         impacts(i)%volume = sum(volume)
         deallocate(price, volume)
      end do
   end function trade_impacts

   pure integer function find_event(events, event_id)
      type(event_t), intent(in) :: events(:)
      integer, intent(in) :: event_id
      integer :: i
      find_event = 0
      do i = 1, size(events)
         if (events(i)%event_id == event_id) then
            find_event = i
            return
         end if
      end do
   end function find_event

   subroutine sort_trades_time(trades)
      type(trade_t), intent(inout) :: trades(:)
      type(trade_t) :: key
      integer :: i, j
      do i = 2, size(trades)
         key = trades(i)
         j = i - 1
         do while (j >= 1)
            if (trades(j)%timestamp_ms <= key%timestamp_ms) exit
            trades(j+1) = trades(j)
            j = j - 1
         end do
         trades(j+1) = key
      end do
   end subroutine sort_trades_time

   subroutine swap_trade_roles(trade, taker_event)
      type(trade_t), intent(inout) :: trade
      type(event_t), intent(in) :: taker_event
      integer :: event_tmp
      integer(i8) :: id_tmp
      trade%price = taker_event%price
      trade%direction = -trade%direction
      event_tmp = trade%maker_event_id
      trade%maker_event_id = trade%taker_event_id
      trade%taker_event_id = event_tmp
      id_tmp = trade%maker_id
      trade%maker_id = trade%taker_id
      trade%taker_id = id_tmp
   end subroutine swap_trade_roles

end module ob_trades
