! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_data
   use pinstimation_kinds, only : dp, i8
   use pinstimation_types, only : trade_counts
   implicit none
   private
   integer, parameter, public :: trade_sell = -1, trade_unresolved = 0, trade_buy = 1
   public :: classify_trades, aggregate_classifications

contains

   subroutine classify_trades(timestamp, price, bid, ask, method, classification, timelag, status)
      real(dp), intent(in) :: timestamp(:), price(:), bid(:), ask(:)
      character(len=*), intent(in) :: method
      integer, allocatable, intent(out) :: classification(:)
      real(dp), intent(in), optional :: timelag
      integer, intent(out), optional :: status
      integer, allocatable :: tick(:)
      integer :: i, j, n
      real(dp) :: lag, threshold, midpoint, difference
      character(len=8) :: selected
      if (present(status)) status = 0
      n = size(timestamp)
      allocate(classification(n), tick(n))
      classification = trade_unresolved
      tick = trade_unresolved
      if (size(price) /= n .or. size(bid) /= n .or. size(ask) /= n) then
         if (present(status)) status = 1
         return
      end if
      do i = 2, n
         if (price(i) > price(i-1)) then
            tick(i) = trade_buy
         else if (price(i) < price(i-1)) then
            tick(i) = trade_sell
         else
            tick(i) = tick(i-1)
         end if
      end do
      selected = uppercase(adjustl(method))
      lag = 0.0_dp
      if (present(timelag)) lag = timelag
      if (trim(selected) == 'TICK') then
         classification = tick
         return
      end if
      do i = 1, n
         threshold = timestamp(i) - lag
         j = lag_index(timestamp, threshold, lag >= 0.0_dp)
         if (j == 0) cycle
         select case (trim(selected))
         case ('EMO')
            if (abs(price(i) - ask(j)) <= 8.0_dp*epsilon(price(i))*max(1.0_dp,abs(price(i)),abs(ask(j)))) then
               classification(i) = trade_buy
            else if (abs(price(i) - bid(j)) <= 8.0_dp*epsilon(price(i))*max(1.0_dp,abs(price(i)),abs(bid(j)))) then
               classification(i) = trade_sell
            else
               classification(i) = tick(i)
            end if
         case ('LR','QUOTE')
            midpoint = 0.5_dp*(bid(j) + ask(j))
            difference = price(i) - midpoint
            if (difference > 0.0_dp) then
               classification(i) = trade_buy
            else if (difference < 0.0_dp) then
               classification(i) = trade_sell
            else if (trim(selected) == 'LR') then
               classification(i) = tick(i)
            end if
         case default
            classification(i) = tick(i)
         end select
      end do
   end subroutine classify_trades

   subroutine aggregate_classifications(group, classification, counts, status)
      integer, intent(in) :: group(:), classification(:)
      type(trade_counts), intent(out) :: counts
      integer, intent(out), optional :: status
      integer :: i, ng
      if (present(status)) status = 0
      if (size(group) /= size(classification) .or. size(group) == 0 .or. minval(group) < 1) then
         allocate(counts%buys(0), counts%sells(0))
         if (present(status)) status = 1
         return
      end if
      ng = maxval(group)
      allocate(counts%buys(ng), counts%sells(ng))
      counts%buys = 0_i8; counts%sells = 0_i8
      do i = 1, size(group)
         if (classification(i) == trade_buy) counts%buys(group(i)) = counts%buys(group(i)) + 1_i8
         if (classification(i) == trade_sell) counts%sells(group(i)) = counts%sells(group(i)) + 1_i8
      end do
   end subroutine aggregate_classifications

   integer function lag_index(timestamp, threshold, backward) result(index)
      real(dp), intent(in) :: timestamp(:), threshold
      logical, intent(in) :: backward
      integer :: lo, hi, mid
      index = 0
      lo = 1; hi = size(timestamp)
      if (backward) then
         if (threshold < timestamp(1)) return
         do while (lo <= hi)
            mid = (lo + hi)/2
            if (timestamp(mid) <= threshold) then
               index = mid; lo = mid + 1
            else
               hi = mid - 1
            end if
         end do
      else
         if (threshold > timestamp(size(timestamp))) return
         index = size(timestamp) + 1
         do while (lo <= hi)
            mid = (lo + hi)/2
            if (timestamp(mid) >= threshold) then
               index = mid; hi = mid - 1
            else
               lo = mid + 1
            end if
         end do
         if (index > size(timestamp)) index = 0
      end if
   end function lag_index

   pure function uppercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) then
            out(i:i) = achar(code - 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function uppercase

end module pinstimation_data
