! SPDX-License-Identifier: GPL-2.0-or-later
program event_matching
   use ob_analytics
   implicit none
   type(event_t)::events(4)
   integer::i
   do i=1,4
      events(i)%event_id=i
      events(i)%id=int(i,i8)
      events(i)%fill=100.0_dp
   end do
   events%timestamp_ms=[0_i8,10_i8,10000_i8,10010_i8]
   events%side=[side_bid,side_ask,side_bid,side_ask]
   call event_match(events,1000_i8)
   do i=1,4
      print '(i0,1x,i0)',events(i)%event_id,events(i)%matching_event
   end do
end program event_matching
