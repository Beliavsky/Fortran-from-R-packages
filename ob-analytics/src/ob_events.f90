! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_events
   use ob_kinds, only : dp, i8
   use ob_types
   use ob_alignment, only : similarity_matrix_time, needleman_wunsch
   implicit none
   private
   public :: event_match, set_order_types, order_aggressiveness

contains

   subroutine event_match(events, cutoff_ms, status, message)
      type(event_t), intent(inout) :: events(:)
      integer(i8), intent(in), optional :: cutoff_ms
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: message
      integer(i8) :: cutoff
      integer, allocatable :: bid_candidates(:), ask_candidates(:)
      integer, allocatable :: pair_bid(:), pair_ask(:)
      real(dp), allocatable :: volumes(:)
      integer :: i, j, k, nvol, nb, na
      integer, allocatable :: bids(:), asks(:), nearest(:), alignment(:,:)
      real(dp), allocatable :: score(:,:)
      integer(i8), allocatable :: bid_time(:), ask_time(:)
      real(dp) :: best_distance, distance
      integer :: best_j
      logical :: duplicated

      cutoff = 5000_i8
      if (present(cutoff_ms)) cutoff = cutoff_ms
      if (present(status)) status = 0
      if (present(message)) message = ''
      events%matching_event = 0
      if (cutoff < 0_i8) then
         call fail(1, 'event_match: cutoff_ms must be nonnegative')
         return
      end if

      bid_candidates = pack([(i, i=1,size(events))], &
         events%side == side_bid .and. events%fill /= 0.0_dp)
      ask_candidates = pack([(i, i=1,size(events))], &
         events%side == side_ask .and. events%fill /= 0.0_dp)
      if (size(bid_candidates) == 0 .or. size(ask_candidates) == 0) return

      allocate(volumes(size(bid_candidates)))
      nvol = 0
      do i = 1, size(bid_candidates)
         if (.not. any(events(ask_candidates)%fill == events(bid_candidates(i))%fill)) cycle
         if (.not. any(volumes(1:nvol) == events(bid_candidates(i))%fill)) then
            nvol = nvol + 1
            volumes(nvol) = events(bid_candidates(i))%fill
         end if
      end do
      if (nvol == 0) return
      if (nvol < size(volumes)) volumes = volumes(1:nvol)
      allocate(pair_bid(0), pair_ask(0))

      do k = 1, nvol
         bids = pack(bid_candidates, events(bid_candidates)%fill == volumes(k))
         asks = pack(ask_candidates, events(ask_candidates)%fill == volumes(k))
         call sort_event_indices_time(events, bids)
         call sort_event_indices_time(events, asks)
         nb = size(bids)
         na = size(asks)
         allocate(nearest(nb))
         nearest = 0
         do i = 1, nb
            best_distance = huge(1.0_dp)
            best_j = 0
            do j = 1, na
               distance = real(abs(events(bids(i))%timestamp_ms-events(asks(j))%timestamp_ms), dp)
               if (distance < best_distance) then
                  best_distance = distance
                  best_j = j
               end if
            end do
            if (best_distance <= real(cutoff,dp)) nearest(i) = best_j
         end do

         duplicated = has_duplicate_nonzero(nearest)
         if (count(nearest == 0) > 1) duplicated = .true.
         if (.not. duplicated) then
            do i = 1, nb
               if (nearest(i) > 0) call append_pair(pair_bid, pair_ask, bids(i), asks(nearest(i)))
            end do
         else
            bid_time = events(bids)%timestamp_ms
            ask_time = events(asks)%timestamp_ms
            score = similarity_matrix_time(bid_time, ask_time, real(cutoff,dp))
            alignment = needleman_wunsch(score, -1.0_dp)
            do i = 1, size(alignment,1)
               if (abs(events(bids(alignment(i,1)))%timestamp_ms - &
                       events(asks(alignment(i,2)))%timestamp_ms) <= cutoff) then
                  call append_pair(pair_bid, pair_ask, bids(alignment(i,1)), asks(alignment(i,2)))
               end if
            end do
         end if
         deallocate(nearest)
      end do

      if (has_duplicate_nonzero(pair_bid) .or. has_duplicate_nonzero(pair_ask)) then
         call fail(2, 'event_match: matching produced duplicate events')
         return
      end if
      do i = 1, size(pair_bid)
         events(pair_bid(i))%matching_event = events(pair_ask(i))%event_id
         events(pair_ask(i))%matching_event = events(pair_bid(i))%event_id
      end do

   contains
      subroutine fail(code, text)
         integer, intent(in) :: code
         character(len=*), intent(in) :: text
         if (present(status)) status = code
         if (present(message)) message = text
      end subroutine fail
   end subroutine event_match

   subroutine set_order_types(events, trades, unidentified)
      type(event_t), intent(inout) :: events(:)
      type(trade_t), intent(in) :: trades(:)
      integer, intent(out), optional :: unidentified
      integer(i8), allocatable :: ids(:), maker_ids(:), taker_ids(:)
      integer :: i, j
      logical :: has_created, has_changed, has_deleted, pacman, flashed
      real(dp) :: created_volume, deleted_volume
      logical :: is_maker, is_taker

      ids = unique_event_ids(events)
      allocate(maker_ids(size(trades)), taker_ids(size(trades)))
      do i = 1, size(trades)
         maker_ids(i) = event_order_id(events, trades(i)%maker_event_id)
         taker_ids(i) = event_order_id(events, trades(i)%taker_event_id)
      end do
      events%order_type = type_unknown

      do i = 1, size(ids)
         has_created = .false.
         has_changed = .false.
         has_deleted = .false.
         pacman = .false.
         created_volume = 0.0_dp
         deleted_volume = 0.0_dp
         do j = 1, size(events)
            if (events(j)%id /= ids(i)) cycle
            select case (events(j)%action)
            case (action_created)
               has_created = .true.
               created_volume = events(j)%volume
            case (action_changed)
               has_changed = .true.
            case (action_deleted)
               has_deleted = .true.
               deleted_volume = events(j)%volume
            end select
         end do
         pacman = order_changes_price(events, ids(i))
         is_maker = any(maker_ids == ids(i))
         is_taker = any(taker_ids == ids(i))
         flashed = has_created .and. has_deleted .and. .not. has_changed .and. &
            deleted_volume == created_volume

         if (pacman) then
            call assign_type(events, ids(i), type_pacman)
         else if (is_taker .and. is_maker) then
            call assign_type(events, ids(i), type_market_limit)
         else if (is_taker) then
            call assign_type(events, ids(i), type_market)
         else if (flashed) then
            call assign_type(events, ids(i), type_flashed_limit)
         else if ((has_created .and. .not. has_changed .and. .not. has_deleted) .or. &
                  (is_maker .and. .not. is_taker)) then
            call assign_type(events, ids(i), type_resting_limit)
         end if
      end do
      if (present(unidentified)) unidentified = count(events%order_type == type_unknown)
   end subroutine set_order_types

   subroutine order_aggressiveness(events, summary, status)
      type(event_t), intent(inout) :: events(:)
      type(depth_summary_t), intent(in) :: summary
      integer, intent(out), optional :: status
      integer, allocatable :: selected(:)
      integer :: side, i, j, previous_summary
      real(dp) :: best, direction
      if (present(status)) status = 0
      events%has_aggressiveness = .false.
      do side = side_bid, side_ask
         selected = pack([(i,i=1,size(events))], events%side == side .and. &
            events%action /= action_changed .and. &
            (events%order_type == type_flashed_limit .or. &
             events%order_type == type_resting_limit))
         call sort_event_indices_time(events, selected)
         direction = merge(1.0_dp, -1.0_dp, side == side_bid)
         do j = 2, size(selected)
            previous_summary = find_timestamp(summary%timestamp_ms, events(selected(j-1))%timestamp_ms)
            if (previous_summary == 0) then
               if (present(status)) status = 1
               cycle
            end if
            if (side == side_bid) then
               best = summary%best_bid_price(previous_summary)
            else
               best = summary%best_ask_price(previous_summary)
            end if
            if (best > 0.0_dp) then
               events(selected(j))%aggressiveness_bps = &
                  10000.0_dp*direction*(events(selected(j))%price-best)/best
               events(selected(j))%has_aggressiveness = .true.
            end if
         end do
      end do
   end subroutine order_aggressiveness

   subroutine sort_event_indices_time(events, idx)
      type(event_t), intent(in) :: events(:)
      integer, intent(inout) :: idx(:)
      integer :: i, j, key
      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (events(idx(j))%timestamp_ms < events(key)%timestamp_ms) exit
            if (events(idx(j))%timestamp_ms == events(key)%timestamp_ms .and. &
                events(idx(j))%event_id <= events(key)%event_id) exit
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end subroutine sort_event_indices_time

   pure logical function has_duplicate_nonzero(v)
      integer, intent(in) :: v(:)
      integer :: i
      has_duplicate_nonzero = .false.
      do i = 2, size(v)
         if (v(i) /= 0 .and. any(v(1:i-1) == v(i))) then
            has_duplicate_nonzero = .true.
            return
         end if
      end do
   end function has_duplicate_nonzero

   subroutine append_pair(a, b, av, bv)
      integer, allocatable, intent(inout) :: a(:), b(:)
      integer, intent(in) :: av, bv
      integer, allocatable :: ta(:), tb(:)
      integer :: n
      n = size(a)
      allocate(ta(n+1), tb(n+1))
      if (n > 0) then
         ta(1:n) = a
         tb(1:n) = b
      end if
      ta(n+1) = av
      tb(n+1) = bv
      call move_alloc(ta, a)
      call move_alloc(tb, b)
   end subroutine append_pair

   function unique_event_ids(events) result(ids)
      type(event_t), intent(in) :: events(:)
      integer(i8), allocatable :: ids(:)
      integer :: i, n
      allocate(ids(size(events)))
      n = 0
      do i = 1, size(events)
         if (.not. any(ids(1:n) == events(i)%id)) then
            n = n + 1
            ids(n) = events(i)%id
         end if
      end do
      if (n < size(ids)) ids = ids(1:n)
   end function unique_event_ids

   pure integer(i8) function event_order_id(events, event_id)
      type(event_t), intent(in) :: events(:)
      integer, intent(in) :: event_id
      integer :: i
      event_order_id = 0_i8
      do i = 1, size(events)
         if (events(i)%event_id == event_id) then
            event_order_id = events(i)%id
            return
         end if
      end do
   end function event_order_id


   pure logical function order_changes_price(events, id)
      type(event_t), intent(in) :: events(:)
      integer(i8), intent(in) :: id
      integer :: i
      real(dp) :: previous_price
      logical :: have_previous
      order_changes_price = .false.
      have_previous = .false.
      previous_price = 0.0_dp
      do i = 1, size(events)
         if (events(i)%id /= id) cycle
         if (have_previous .and. events(i)%price /= previous_price) then
            order_changes_price = .true.
            return
         end if
         previous_price = events(i)%price
         have_previous = .true.
      end do
   end function order_changes_price

   subroutine assign_type(events, id, order_type)
      type(event_t), intent(inout) :: events(:)
      integer(i8), intent(in) :: id
      integer, intent(in) :: order_type
      integer :: i
      do i = 1, size(events)
         if (events(i)%id == id) events(i)%order_type = order_type
      end do
   end subroutine assign_type

   pure integer function find_timestamp(timestamps, target)
      integer(i8), intent(in) :: timestamps(:), target
      integer :: i
      find_timestamp = 0
      do i = 1, size(timestamps)
         if (timestamps(i) == target) then
            find_timestamp = i
            return
         end if
      end do
   end function find_timestamp

end module ob_events
