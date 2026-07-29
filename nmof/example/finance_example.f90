! SPDX-License-Identifier: GPL-3.0-only
program finance_example
   use nmof
   implicit none
   real(dp) :: cashflows(6),times(6),all_cashflows(7),all_times(7),price,yield_value
   type(option_result) :: call
   integer :: status

   cashflows=[5.0_dp,5.0_dp,5.0_dp,5.0_dp,5.0_dp,105.0_dp]
   times=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   price=vanilla_bond(cashflows,times,yields=[0.03_dp,0.03_dp,0.03_dp,0.03_dp,0.03_dp,0.03_dp])
   all_cashflows=[-price,cashflows]
   all_times=[0.0_dp,times]
   yield_value=yield_to_maturity(all_cashflows,all_times,tol=1.0e-10_dp,status=status)
   call=vanilla_option_european(100.0_dp,105.0_dp,0.5_dp,0.03_dp,0.01_dp,0.25_dp**2,'call')

   write(*,'(a,f12.6)') 'Bond price: ',price
   write(*,'(a,f12.8)') 'Recovered yield: ',yield_value
   write(*,'(a,f12.6)') 'Option value: ',call%value
end program finance_example
