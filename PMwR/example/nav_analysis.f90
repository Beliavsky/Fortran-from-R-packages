program nav_analysis
   use pmwr, only : dp, nav_summary_result, drawdown_table, summarize_nav, compute_drawdowns
   implicit none
   type(nav_summary_result) :: summary
   type(drawdown_table) :: dd
   real(dp) :: nav(8)
   integer :: i

   nav = [100.0_dp, 106.0_dp, 103.0_dp, 111.0_dp, 95.0_dp, 99.0_dp, 108.0_dp, 116.0_dp]
   call summarize_nav(nav, summary, periods_per_year=12.0_dp)
   call compute_drawdowns(nav, dd)

   print '(a,f8.3,a)', 'Total return: ', 100.0_dp*summary%total_return, '%'
   print '(a,f8.3,a)', 'Maximum drawdown: ', 100.0_dp*summary%max_drawdown, '%'
   do i = 1, size(dd%depth)
      print '(a,i0,a,i0,a,i0,a,f8.3,a)', 'Peak ', dd%peak(i), ', trough ', dd%trough(i), &
            ', recovery ', dd%recover(i), ', depth ', 100.0_dp*dd%depth(i), '%'
   end do
end program nav_analysis
