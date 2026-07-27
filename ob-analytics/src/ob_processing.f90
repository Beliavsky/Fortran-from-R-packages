! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_processing
   use ob_kinds, only : dp, i8
   use ob_types
   use ob_io, only : read_event_csv
   use ob_events, only : event_match, set_order_types, order_aggressiveness
   use ob_trades, only : match_trades
   use ob_depth, only : price_level_volume, depth_metrics
   implicit none
   private
   public :: process_data

contains

   subroutine process_data(csv_file, result, price_digits, volume_digits, &
      cutoff_ms, warmup_ms, remove_zombies, status, message)
      character(len=*), intent(in) :: csv_file
      type(processing_result_t), intent(out) :: result
      integer, intent(in), optional :: price_digits, volume_digits
      integer(i8), intent(in), optional :: cutoff_ms, warmup_ms
      logical, intent(in), optional :: remove_zombies
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: message
      integer :: io_status, match_status, pd, vd
      integer(i8) :: cutoff, warmup, offset
      logical :: drop_zombies
      character(len=:), allocatable :: io_message, match_message
      integer(i8), allocatable :: zombie_ids(:)
      logical, allocatable :: keep(:)

      pd=2; vd=8; cutoff=5000_i8; warmup=60000_i8; drop_zombies=.true.
      if(present(price_digits)) pd=price_digits
      if(present(volume_digits)) vd=volume_digits
      if(present(cutoff_ms)) cutoff=cutoff_ms
      if(present(warmup_ms)) warmup=warmup_ms
      if(present(remove_zombies)) drop_zombies=remove_zombies
      if(present(status)) status=0
      if(present(message)) message=''

      call read_event_csv(csv_file,result%events,pd,vd,io_status,io_message)
      if(io_status/=0) then
         call fail(io_status,io_message); return
      end if
      call event_match(result%events,cutoff,match_status,match_message)
      if(match_status/=0) then
         call fail(100+match_status,match_message); return
      end if
      result%trades=match_trades(result%events)
      call set_order_types(result%events,result%trades)
      if(drop_zombies) then
         zombie_ids=find_zombie_ids(result%events,result%trades)
         if(size(zombie_ids)>0) then
            allocate(keep(size(result%events)))
            keep=.true.
            keep=.true.
            call mark_non_zombies(result%events,zombie_ids,keep)
            result%events=pack(result%events,keep)
         end if
      end if
      result%depth=price_level_volume(result%events)
      result%depth_summary=depth_metrics(result%depth)
      call order_aggressiveness(result%events,result%depth_summary)
      if(size(result%events)>0 .and. warmup>0_i8) then
         offset=minval(result%events%timestamp_ms)+warmup
         call trim_summary(result%depth_summary,offset)
      end if

   contains
      subroutine fail(code,text)
         integer,intent(in)::code
         character(len=*),intent(in)::text
         if(present(status)) status=code
         if(present(message)) message=text
      end subroutine fail
   end subroutine process_data

   function find_zombie_ids(events,trades) result(ids)
      type(event_t),intent(in)::events(:)
      type(trade_t),intent(in)::trades(:)
      integer(i8),allocatable::ids(:)
      integer(i8),allocatable::candidate(:)
      integer::i,j,n,last
      logical::deleted,zombie
      allocate(candidate(size(events)),ids(size(events)))
      n=0
      do i=1,size(events)
         if(any(candidate(1:n)==events(i)%id)) cycle
         deleted=any(events%id==events(i)%id .and. events%action==action_deleted)
         if(.not.deleted) then
            n=n+1; candidate(n)=events(i)%id
         end if
      end do
      ids=0_i8; n=0
      do i=1,count(candidate/=0_i8)
         last=0
         do j=1,size(events)
            if(events(j)%id==candidate(i)) then
               if(last==0) then
                  last=j
               else if(events(j)%timestamp_ms>=events(last)%timestamp_ms) then
                  last=j
               end if
            end if
         end do
         zombie=.false.
         if(last>0) then
            do j=1,size(trades)
               if(trades(j)%timestamp_ms<events(last)%timestamp_ms) cycle
               if(events(last)%side==side_bid .and. trades(j)%direction==trade_sell .and. &
                  trades(j)%price<events(last)%price) zombie=.true.
               if(events(last)%side==side_ask .and. trades(j)%direction==trade_buy .and. &
                  trades(j)%price>events(last)%price) zombie=.true.
            end do
         end if
         if(zombie) then
            n=n+1; ids(n)=candidate(i)
         end if
      end do
      if(n<size(ids)) ids=ids(1:n)
   end function find_zombie_ids

   subroutine mark_non_zombies(events,zombie_ids,keep)
      type(event_t),intent(in)::events(:)
      integer(i8),intent(in)::zombie_ids(:)
      logical,intent(out)::keep(:)
      integer::i
      do i=1,size(events)
         keep(i)=.not.any(zombie_ids==events(i)%id)
      end do
   end subroutine mark_non_zombies

   subroutine trim_summary(summary,offset)
      type(depth_summary_t),intent(inout)::summary
      integer(i8),intent(in)::offset
      logical,allocatable::keep(:)
      allocate(keep(size(summary%timestamp_ms)))
      keep=summary%timestamp_ms>=offset
      summary%timestamp_ms=pack(summary%timestamp_ms,keep)
      summary%best_bid_price=pack(summary%best_bid_price,keep)
      summary%best_bid_volume=pack(summary%best_bid_volume,keep)
      summary%best_ask_price=pack(summary%best_ask_price,keep)
      summary%best_ask_volume=pack(summary%best_ask_volume,keep)
      summary%bid_volume=pack_rows(summary%bid_volume,keep)
      summary%ask_volume=pack_rows(summary%ask_volume,keep)
   end subroutine trim_summary

   function pack_rows(x,keep) result(y)
      real(dp),intent(in)::x(:,:)
      logical,intent(in)::keep(:)
      real(dp),allocatable::y(:,:)
      integer::i,j
      allocate(y(count(keep),size(x,2)))
      j=0
      do i=1,size(keep)
         if(keep(i)) then
            j=j+1; y(j,:)=x(i,:)
         end if
      end do
   end function pack_rows

end module ob_processing
