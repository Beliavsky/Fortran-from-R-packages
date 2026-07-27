! SPDX-License-Identifier: GPL-2.0-or-later
program ob_analytics_demo
   use ob_analytics
   implicit none
   type(processing_result_t)::result
   integer::status
   character(len=:),allocatable::message
   call process_data('example/data/orders.csv',result,warmup_ms=0_i8, &
      remove_zombies=.false.,status=status,message=message)
   if(status/=0) then
      print '(a)',message
      error stop 1
   end if
   print '(a,i0)', 'events: ',size(result%events)
   print '(a,i0)', 'trades: ',size(result%trades)
   print '(a,i0)', 'depth updates: ',size(result%depth)
   if(size(result%trades)>0) then
      print '(a,f8.2,a,f8.2)', 'first trade price: ',result%trades(1)%price, &
         ' volume: ',result%trades(1)%volume
   end if
end program ob_analytics_demo
