! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_depth
   use ob_kinds, only : dp, i8
   use ob_types
   implicit none
   private
   public :: price_level_volume, filter_depth, depth_metrics, get_spread

contains

   function price_level_volume(events) result(depth)
      type(event_t), intent(in) :: events(:)
      type(depth_update_t), allocatable :: depth(:)
      type(depth_update_t), allocatable :: bid(:), ask(:)
      bid = directional_depth(events, side_bid)
      ask = directional_depth(events, side_ask)
      allocate(depth(size(bid)+size(ask)))
      if (size(bid) > 0) depth(1:size(bid)) = bid
      if (size(ask) > 0) depth(size(bid)+1:) = ask
      call sort_depth_time(depth)
   end function price_level_volume

   function directional_depth(events, side) result(depth)
      type(event_t), intent(in) :: events(:)
      integer, intent(in) :: side
      type(depth_update_t), allocatable :: depth(:)
      type(depth_update_t), allocatable :: deltas(:)
      integer(i8), allocatable :: added_ids(:)
      integer :: i, n
      real(dp) :: delta
      logical :: include

      allocate(added_ids(size(events)), deltas(size(events)*2))
      n = 0
      added_ids = 0_i8
      do i = 1, size(events)
         if (events(i)%side /= side) cycle
         if (events(i)%order_type == type_pacman .or. events(i)%order_type == type_market) cycle
         if (events(i)%action == action_created .or. &
             (events(i)%action == action_changed .and. events(i)%fill == 0.0_dp)) then
            if (.not. any(added_ids == events(i)%id)) added_ids(count(added_ids /= 0_i8)+1) = events(i)%id
         end if
      end do

      do i = 1, size(events)
         if (events(i)%side /= side) cycle
         if (events(i)%order_type == type_pacman .or. events(i)%order_type == type_market) cycle
         include = .false.
         delta = 0.0_dp
         if (events(i)%action == action_created .or. &
             (events(i)%action == action_changed .and. events(i)%fill == 0.0_dp)) then
            include = .true.
            delta = events(i)%volume
         else if (events(i)%action == action_deleted .and. events(i)%volume > 0.0_dp .and. &
                  any(added_ids == events(i)%id)) then
            include = .true.
            delta = -events(i)%volume
         end if
         if (events(i)%fill > 0.0_dp .and. any(added_ids == events(i)%id)) then
            if (include) then
               call append_delta(deltas, n, events(i), delta)
            end if
            include = .true.
            delta = -events(i)%fill
         end if
         if (include) call append_delta(deltas, n, events(i), delta)
      end do
      if (n == 0) then
         allocate(depth(0))
         return
      end if
      deltas = deltas(1:n)
      call sort_depth_price_time(deltas)
      do i = 1, n
         if (i > 1) then
            if (deltas(i)%price == deltas(i-1)%price) then
               deltas(i)%volume = max(0.0_dp, deltas(i-1)%volume + deltas(i)%volume)
            else
               deltas(i)%volume = max(0.0_dp, deltas(i)%volume)
            end if
         else
            deltas(i)%volume = max(0.0_dp, deltas(i)%volume)
         end if
      end do
      call move_alloc(deltas, depth)
   end function directional_depth

   function filter_depth(depth, from_ms, to_ms) result(filtered)
      type(depth_update_t), intent(in) :: depth(:)
      integer(i8), intent(in) :: from_ms, to_ms
      type(depth_update_t), allocatable :: filtered(:)
      type(depth_update_t), allocatable :: work(:)
      integer :: i, j, n, latest
      logical :: already
      if (to_ms <= from_ms) then
         allocate(filtered(0))
         return
      end if
      allocate(work(2*size(depth)+1))
      n = 0
      do i = 1, size(depth)
         if (depth(i)%timestamp_ms > from_ms) cycle
         latest = i
         do j = i+1, size(depth)
            if (same_level(depth(j), depth(i)) .and. depth(j)%timestamp_ms <= from_ms .and. &
                depth(j)%timestamp_ms >= depth(latest)%timestamp_ms) latest = j
         end do
         already = .false.
         do j = 1, i-1
            if (same_level(depth(j),depth(i))) then
               already = .true.
               exit
            end if
         end do
         if (.not. already .and. depth(latest)%volume > 0.0_dp) then
            n = n + 1
            work(n) = depth(latest)
            work(n)%timestamp_ms = from_ms
         end if
      end do
      do i = 1, size(depth)
         if (depth(i)%timestamp_ms > from_ms .and. depth(i)%timestamp_ms < to_ms) then
            n = n + 1
            work(n) = depth(i)
         end if
      end do
      do i = 1, n
         already = .false.
         do j = i+1, n
            if (same_level(work(j),work(i))) then
               already = .true.
               exit
            end if
         end do
         if (.not. already .and. work(i)%volume > 0.0_dp) then
            n = n + 1
            work(n) = work(i)
            work(n)%timestamp_ms = to_ms
            work(n)%volume = 0.0_dp
         end if
      end do
      if (n == 0) then
         allocate(filtered(0))
      else
         work = work(1:n)
         call sort_depth_price_time(work)
         call move_alloc(work, filtered)
      end if
   end function filter_depth

   function depth_metrics(depth, bps, bins) result(summary)
      type(depth_update_t), intent(in) :: depth(:)
      integer, intent(in), optional :: bps, bins
      type(depth_summary_t) :: summary
      integer :: bps_value, bins_value, i
      integer, allocatable :: bid_cents(:), ask_cents(:)
      real(dp), allocatable :: bid_vol(:), ask_vol(:)
      integer :: cents, best_bid_cents, best_ask_cents

      bps_value = 25
      bins_value = 20
      if (present(bps)) bps_value = bps
      if (present(bins)) bins_value = bins
      if (bps_value <= 0 .or. bins_value <= 0) error stop 'depth_metrics: bps and bins must be positive'
      summary%bps = bps_value
      summary%bins = bins_value
      allocate(summary%timestamp_ms(size(depth)), summary%best_bid_price(size(depth)), &
         summary%best_bid_volume(size(depth)), summary%bid_volume(size(depth),bins_value), &
         summary%best_ask_price(size(depth)), summary%best_ask_volume(size(depth)), &
         summary%ask_volume(size(depth),bins_value))
      summary%best_bid_price = 0.0_dp
      summary%best_bid_volume = 0.0_dp
      summary%bid_volume = 0.0_dp
      summary%best_ask_price = 0.0_dp
      summary%best_ask_volume = 0.0_dp
      summary%ask_volume = 0.0_dp
      allocate(bid_cents(0), ask_cents(0), bid_vol(0), ask_vol(0))

      do i = 1, size(depth)
         summary%timestamp_ms(i) = depth(i)%timestamp_ms
         cents = nint(100.0_dp*depth(i)%price)
         best_bid_cents = best_level(bid_cents,bid_vol,.true.)
         best_ask_cents = best_level(ask_cents,ask_vol,.false.)
         if (depth(i)%side == side_ask) then
            if (best_bid_cents == 0 .or. cents > best_bid_cents) &
               call set_level(ask_cents, ask_vol, cents, depth(i)%volume)
         else
            if (best_ask_cents == 0 .or. cents < best_ask_cents) &
               call set_level(bid_cents, bid_vol, cents, depth(i)%volume)
         end if
         best_bid_cents = best_level(bid_cents,bid_vol,.true.)
         best_ask_cents = best_level(ask_cents,ask_vol,.false.)
         if (best_bid_cents > 0) then
            summary%best_bid_price(i) = real(best_bid_cents,dp)/100.0_dp
            summary%best_bid_volume(i) = level_volume(bid_cents,bid_vol,best_bid_cents)
            summary%bid_volume(i,:) = volume_bins(bid_cents,bid_vol,best_bid_cents,bps_value,bins_value,.true.)
         end if
         if (best_ask_cents > 0) then
            summary%best_ask_price(i) = real(best_ask_cents,dp)/100.0_dp
            summary%best_ask_volume(i) = level_volume(ask_cents,ask_vol,best_ask_cents)
            summary%ask_volume(i,:) = volume_bins(ask_cents,ask_vol,best_ask_cents,bps_value,bins_value,.false.)
         end if
      end do
   end function depth_metrics

   function get_spread(summary) result(spread)
      type(depth_summary_t), intent(in) :: summary
      type(spread_t) :: spread
      logical, allocatable :: keep(:)
      integer :: i, n
      allocate(keep(size(summary%timestamp_ms)))
      keep = .false.
      if (size(keep) > 0) keep(1) = .true.
      do i = 2, size(keep)
         keep(i) = summary%best_bid_price(i) /= summary%best_bid_price(i-1) .or. &
            summary%best_bid_volume(i) /= summary%best_bid_volume(i-1) .or. &
            summary%best_ask_price(i) /= summary%best_ask_price(i-1) .or. &
            summary%best_ask_volume(i) /= summary%best_ask_volume(i-1)
      end do
      n = count(keep)
      allocate(spread%timestamp_ms(n), spread%bid_price(n), spread%bid_volume(n), &
         spread%ask_price(n), spread%ask_volume(n))
      spread%timestamp_ms = pack(summary%timestamp_ms,keep)
      spread%bid_price = pack(summary%best_bid_price,keep)
      spread%bid_volume = pack(summary%best_bid_volume,keep)
      spread%ask_price = pack(summary%best_ask_price,keep)
      spread%ask_volume = pack(summary%best_ask_volume,keep)
   end function get_spread

   subroutine append_delta(deltas, n, event, volume)
      type(depth_update_t), intent(inout) :: deltas(:)
      integer, intent(inout) :: n
      type(event_t), intent(in) :: event
      real(dp), intent(in) :: volume
      n = n + 1
      deltas(n)%timestamp_ms = event%timestamp_ms
      deltas(n)%price = event%price
      deltas(n)%volume = volume
      deltas(n)%side = event%side
   end subroutine append_delta

   subroutine sort_depth_price_time(depth)
      type(depth_update_t), intent(inout) :: depth(:)
      type(depth_update_t) :: key
      integer :: i, j
      do i = 2, size(depth)
         key = depth(i)
         j = i-1
         do while (j >= 1)
            if (depth(j)%price < key%price) exit
            if (depth(j)%price == key%price .and. depth(j)%timestamp_ms <= key%timestamp_ms) exit
            depth(j+1) = depth(j)
            j = j-1
         end do
         depth(j+1) = key
      end do
   end subroutine sort_depth_price_time

   subroutine sort_depth_time(depth)
      type(depth_update_t), intent(inout) :: depth(:)
      type(depth_update_t) :: key
      integer :: i, j
      do i = 2, size(depth)
         key = depth(i)
         j = i-1
         do while (j >= 1)
            if (depth(j)%timestamp_ms <= key%timestamp_ms) exit
            depth(j+1) = depth(j)
            j = j-1
         end do
         depth(j+1) = key
      end do
   end subroutine sort_depth_time

   pure logical function same_level(a,b)
      type(depth_update_t), intent(in) :: a,b
      same_level = a%price == b%price .and. a%side == b%side
   end function same_level

   subroutine set_level(cents, volume, price_cents, value)
      integer, allocatable, intent(inout) :: cents(:)
      real(dp), allocatable, intent(inout) :: volume(:)
      integer, intent(in) :: price_cents
      real(dp), intent(in) :: value
      integer :: i, n
      integer, allocatable :: new_cents(:)
      real(dp), allocatable :: new_volume(:)
      do i = 1, size(cents)
         if (cents(i) == price_cents) then
            volume(i) = value
            return
         end if
      end do
      n = size(cents)
      allocate(new_cents(n+1),new_volume(n+1))
      if (n > 0) then
         new_cents(1:n) = cents
         new_volume(1:n) = volume
      end if
      new_cents(n+1) = price_cents
      new_volume(n+1) = value
      call move_alloc(new_cents,cents)
      call move_alloc(new_volume,volume)
   end subroutine set_level

   pure integer function best_level(cents, volume, bid)
      integer, intent(in) :: cents(:)
      real(dp), intent(in) :: volume(:)
      logical, intent(in) :: bid
      integer :: i
      best_level = 0
      do i = 1, size(cents)
         if (volume(i) <= 0.0_dp) cycle
         if (best_level == 0) then
            best_level = cents(i)
         else if (bid .and. cents(i) > best_level) then
            best_level = cents(i)
         else if (.not. bid .and. cents(i) < best_level) then
            best_level = cents(i)
         end if
      end do
   end function best_level

   pure real(dp) function level_volume(cents, volume, price_cents)
      integer, intent(in) :: cents(:), price_cents
      real(dp), intent(in) :: volume(:)
      integer :: i
      level_volume = 0.0_dp
      do i = 1, size(cents)
         if (cents(i) == price_cents) then
            level_volume = volume(i)
            return
         end if
      end do
   end function level_volume

   pure function volume_bins(cents, volume, best, bps, bins, bid) result(result)
      integer, intent(in) :: cents(:), best, bps, bins
      real(dp), intent(in) :: volume(:)
      logical, intent(in) :: bid
      real(dp) :: result(bins)
      integer :: boundary, nlevels, pos, bin, i, end_level
      result = 0.0_dp
      if (bid) then
         end_level = nint((1.0_dp-real(bps*bins,dp)*0.0001_dp)*real(best,dp))
         nlevels = best-end_level+1
      else
         end_level = nint((1.0_dp+real(bps*bins,dp)*0.0001_dp)*real(best,dp))
         nlevels = end_level-best+1
      end if
      nlevels = max(nlevels,1)
      do i = 1, size(cents)
         if (volume(i) <= 0.0_dp) cycle
         if (bid) then
            if (cents(i) > best .or. cents(i) < end_level) cycle
            pos = best-cents(i)+1
         else
            if (cents(i) < best .or. cents(i) > end_level) cycle
            pos = cents(i)-best+1
         end if
         do bin = 1, bins
            boundary = ceiling(real(bin*nlevels,dp)/real(bins,dp))
            if (pos <= boundary) then
               result(bin) = result(bin)+volume(i)
               exit
            end if
         end do
      end do
   end function volume_bins

end module ob_depth
