program demo_fcl
   use fcl
   implicit none
   type(fixed_bond_type) :: bond
   type(bond_value_type) :: value
   type(return_series_type) :: returns
   real(dp), allocatable :: cumulative(:), dietz_values(:)

   bond%value_date = make_date(2021, 2, 1)
   bond%maturity_date = make_date(2030, 2, 1)
   bond%redemption_value = 100.0_dp
   bond%coupon_rate = 0.03_dp
   bond%coupon_frequency = 1
   value = bond%value(make_date(2022, 2, 1), 100.0_dp)

   print '(a,f12.8)', 'YTM:               ', value%ytm
   print '(a,f12.8)', 'Macaulay duration: ', value%macaulay_duration
   print '(a,f12.8)', 'Modified duration: ', value%modified_duration

   call returns%initialize([0, 4, 9], [100.0_dp, 103.0_dp, 110.0_dp], &
      [0.0_dp, 3.0_dp, 7.0_dp])
   cumulative = returns%twrr_cumulative(1, 9)
   dietz_values = returns%dietz(1, 9)
   print '(a,f12.8)', 'Ending TWRR:        ', cumulative(size(cumulative))
   print '(a,f12.8)', 'Ending Dietz:       ', dietz_values(size(dietz_values))
end program demo_fcl
