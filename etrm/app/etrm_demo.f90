! SPDX-License-Identifier: MIT
program etrm_demo
   use etrm
   implicit none

   real(dp) :: futures(8)
   type(strategy_result) :: hedge
   type(strategy_summary) :: stats
   integer :: status, i
   character(len=:), allocatable :: message

   futures = [100.0_dp, 102.0_dp, 105.0_dp, 103.0_dp, 108.0_dp, 112.0_dp, 109.0_dp, 115.0_dp]
   call cppi(10.0_dp, futures, 0.10_dp, 0.12_dp, hedge, status, message, &
      transaction_cost=0.05_dp, integer_trades=.true.)
   if (status /= etrm_ok) error stop message
   call summarize_strategy(hedge, stats)

   print '(a)', "CPPI buyer example"
   print '(a,f8.3)', "Churn rate: ", stats%churn
   print '(a)', " day  market  position  exposed  portfolio"
   do i = 1, size(futures)
      print '(i4,4f10.3)', i, hedge%market(i), hedge%position(i), &
         hedge%exposed(i), hedge%portfolio(i)
   end do
end program etrm_demo
