! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_book
   use ob_kinds, only : dp, i8
   use ob_types
   implicit none
   private
   public :: reconstruct_order_book

contains

   function reconstruct_order_book(events, timestamp_ms, max_levels, bps_range, &
      min_bid, max_ask) result(book)
      type(event_t), intent(in) :: events(:)
      integer(i8), intent(in) :: timestamp_ms
      integer, intent(in), optional :: max_levels
      real(dp), intent(in), optional :: bps_range, min_bid, max_ask
      type(order_book_t) :: book
      integer(i8), allocatable :: ids(:)
      type(event_t), allocatable :: active(:), bids(:), asks(:)
      integer :: i, n, idx, levels
      real(dp) :: pct, min_bid_value, max_ask_value, best

      book%timestamp_ms = timestamp_ms
      pct = 0.0_dp
      min_bid_value = 0.0_dp
      max_ask_value = huge(1.0_dp)
      if (present(bps_range)) pct = bps_range*0.0001_dp
      if (present(min_bid)) min_bid_value = min_bid
      if (present(max_ask)) max_ask_value = max_ask
      ids = unique_ids(events)
      allocate(active(size(ids)))
      n = 0
      do i = 1, size(ids)
         if (.not. has_action_before(events,ids(i),action_created,timestamp_ms)) cycle
         if (has_action_before(events,ids(i),action_deleted,timestamp_ms)) cycle
         idx = latest_active_event(events,ids(i),timestamp_ms)
         if (idx == 0) cycle
         if (events(idx)%order_type == type_market) cycle
         n = n+1
         active(n) = events(idx)
      end do
      if (n < size(active)) active = active(1:n)
      bids = pack(active, active%side == side_bid)
      asks = pack(active, active%side == side_ask)
      call sort_active(bids,.true.)
      call sort_active(asks,.false.)

      if (size(bids) > 0) then
         best = bids(1)%price
         if (pct > 0.0_dp) min_bid_value = best*(1.0_dp-pct)
         bids = pack(bids,bids%price >= min_bid_value)
      end if
      if (size(asks) > 0) then
         best = asks(1)%price
         if (pct > 0.0_dp) max_ask_value = best*(1.0_dp+pct)
         asks = pack(asks,asks%price <= max_ask_value)
      end if
      if (present(max_levels)) then
         levels = max(max_levels,0)
         if (size(bids) > levels) bids = bids(1:levels)
         if (size(asks) > levels) asks = asks(1:levels)
      end if
      book%bids = make_levels(bids,.true.)
      book%asks = make_levels(asks,.false.)
      call reverse_levels(book%asks)
   end function reconstruct_order_book

   function make_levels(active,bid) result(levels)
      type(event_t), intent(in) :: active(:)
      logical, intent(in) :: bid
      type(order_level_t), allocatable :: levels(:)
      real(dp) :: cumulative, first_price
      integer :: i
      allocate(levels(size(active)))
      if (size(active) == 0) return
      first_price = active(1)%price
      cumulative = 0.0_dp
      do i = 1, size(active)
         cumulative = cumulative+active(i)%volume
         levels(i)%id = active(i)%id
         levels(i)%timestamp_ms = active(i)%timestamp_ms
         levels(i)%exchange_timestamp_ms = active(i)%exchange_timestamp_ms
         levels(i)%price = active(i)%price
         levels(i)%volume = active(i)%volume
         levels(i)%liquidity = cumulative
         if (bid) then
            levels(i)%bps = (first_price-active(i)%price)/first_price*10000.0_dp
         else
            levels(i)%bps = (active(i)%price-first_price)/first_price*10000.0_dp
         end if
      end do
   end function make_levels

   subroutine reverse_levels(levels)
      type(order_level_t), intent(inout) :: levels(:)
      type(order_level_t) :: tmp
      integer :: i,n
      n=size(levels)
      do i=1,n/2
         tmp=levels(i); levels(i)=levels(n-i+1); levels(n-i+1)=tmp
      end do
   end subroutine reverse_levels

   subroutine sort_active(active,bid)
      type(event_t), intent(inout) :: active(:)
      logical, intent(in) :: bid
      type(event_t) :: key
      integer :: i,j
      do i=2,size(active)
         key=active(i); j=i-1
         do while(j>=1)
            if (bid) then
               if (active(j)%price > key%price) exit
               if (active(j)%price == key%price .and. active(j)%id <= key%id) exit
            else
               if (active(j)%price < key%price) exit
               if (active(j)%price == key%price .and. active(j)%id <= key%id) exit
            end if
            active(j+1)=active(j); j=j-1
         end do
         active(j+1)=key
      end do
   end subroutine sort_active

   function unique_ids(events) result(ids)
      type(event_t),intent(in)::events(:)
      integer(i8),allocatable::ids(:)
      integer::i,n
      allocate(ids(size(events))); n=0
      do i=1,size(events)
         if(.not.any(ids(1:n)==events(i)%id)) then
            n=n+1; ids(n)=events(i)%id
         end if
      end do
      if(n<size(ids)) ids=ids(1:n)
   end function unique_ids

   pure logical function has_action_before(events,id,action,timestamp_ms)
      type(event_t),intent(in)::events(:)
      integer(i8),intent(in)::id,timestamp_ms
      integer,intent(in)::action
      integer::i
      has_action_before=.false.
      do i=1,size(events)
         if(events(i)%id==id .and. events(i)%action==action .and. &
            events(i)%timestamp_ms<=timestamp_ms) then
            has_action_before=.true.; return
         end if
      end do
   end function has_action_before

   pure integer function latest_active_event(events,id,timestamp_ms)
      type(event_t),intent(in)::events(:)
      integer(i8),intent(in)::id,timestamp_ms
      integer::i
      latest_active_event=0
      do i=1,size(events)
         if(events(i)%id/=id .or. events(i)%timestamp_ms>timestamp_ms) cycle
         if(events(i)%action==action_deleted) cycle
         if(latest_active_event==0) then
            latest_active_event=i
         else if(events(i)%timestamp_ms >= events(latest_active_event)%timestamp_ms) then
            latest_active_event=i
         end if
      end do
   end function latest_active_event

end module ob_book
