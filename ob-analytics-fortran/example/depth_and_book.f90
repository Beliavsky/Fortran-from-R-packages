! SPDX-License-Identifier: GPL-2.0-or-later
program depth_and_book
   use ob_analytics
   implicit none
   type(event_t)::events(3)
   type(depth_update_t),allocatable::depth(:)
   type(depth_summary_t)::summary
   type(order_book_t)::book
   events(1)=event_t(event_id=1,id=1_i8,timestamp_ms=1000_i8, &
      exchange_timestamp_ms=900_i8,price=99.0_dp,volume=10.0_dp, &
      action=action_created,side=side_bid,order_type=type_resting_limit)
   events(2)=event_t(event_id=2,id=2_i8,timestamp_ms=1100_i8, &
      exchange_timestamp_ms=1000_i8,price=101.0_dp,volume=12.0_dp, &
      action=action_created,side=side_ask,order_type=type_resting_limit)
   events(3)=event_t(event_id=3,id=3_i8,timestamp_ms=1200_i8, &
      exchange_timestamp_ms=1100_i8,price=98.0_dp,volume=5.0_dp, &
      action=action_created,side=side_bid,order_type=type_resting_limit)
   depth=price_level_volume(events)
   summary=depth_metrics(depth,bps=50,bins=4)
   book=reconstruct_order_book(events,1500_i8)
   print '(a,f8.2,a,f8.2)', 'best bid: ',summary%best_bid_price(size(depth)), &
      ' best ask: ',summary%best_ask_price(size(depth))
   print '(a,i0,a,i0)', 'book bids: ',size(book%bids),' asks: ',size(book%asks)
end program depth_and_book
