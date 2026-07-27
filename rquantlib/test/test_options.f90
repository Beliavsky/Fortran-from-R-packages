! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program test_options
  use rq_kinds, only: dp
  use rq_options
  use test_support
  implicit none
  type(option_result)::res,call_res,bin_res,asian_res,in_res,out_res
  real(dp)::iv,mc,se,values(3,2),spot,strike,q,r,t,sigma
  integer::status

  res=european_option('call',60.0_dp,65.0_dp,0.0_dp,0.08_dp,0.25_dp,0.3_dp)
  call assert_close(res%value,2.1334_dp,2.0e-4_dp,'European call reference')
  res=european_option('put',100.0_dp,95.0_dp,0.05_dp,0.10_dp,0.5_dp,0.2_dp)
  call assert_close(res%value,2.4648_dp,2.0e-4_dp,'European put reference')
  res=european_option('call',105.0_dp,100.0_dp,0.1_dp,0.1_dp,0.5_dp,0.36_dp)
  call assert_close(res%delta,0.5946_dp,2.0e-4_dp,'European call delta')
  res=european_option('put',105.0_dp,100.0_dp,0.1_dp,0.1_dp,0.5_dp,0.36_dp)
  call assert_close(res%delta,-0.3566_dp,2.0e-4_dp,'European put delta')
  res=european_option('call',55.0_dp,60.0_dp,0.0_dp,0.1_dp,0.75_dp,0.30_dp)
  call assert_close(res%gamma,0.0278_dp,2.0e-4_dp,'European gamma')
  call assert_close(res%vega,18.9358_dp,2.0e-3_dp,'European vega')


  call_res=european_option('call',40.0_dp,40.0_dp,0.0_dp,0.09_dp,0.5_dp,0.3_dp)
  res=european_option('call',40.0_dp,40.0_dp,0.0_dp,0.09_dp,0.5_dp,0.3_dp, &
    [2.0_dp/12.0_dp,5.0_dp/12.0_dp],[0.5_dp,0.5_dp])
  call assert_true(res%value<call_res%value,'Discrete dividends reduce call value')
  bin_res=binary_option('asset','call','american',100.0_dp,105.0_dp,0.01_dp, &
    0.04_dp,0.5_dp,0.25_dp,steps=300)
  call assert_true(bin_res%value>=0.0_dp,'American asset digital')

  res=american_option('call',110.0_dp,100.0_dp,0.1_dp,0.1_dp,0.1_dp,0.15_dp,800)
  call assert_close(res%value,10.0089_dp,3.0e-2_dp,'American call lattice')

  bin_res=binary_option('cash','put','european',100.0_dp,80.0_dp, &
    0.06_dp,0.06_dp,0.75_dp,0.35_dp,10.0_dp)
  call assert_close(bin_res%value,2.6710_dp,3.0e-4_dp,'Cash digital put')

  asian_res=geometric_asian_option('put',80.0_dp,85.0_dp,-0.03_dp, &
    0.05_dp,0.25_dp,0.2_dp)
  call assert_close(asian_res%value,4.6922_dp,4.0e-4_dp,'Geometric Asian put')

  spot=100.0_dp; strike=90.0_dp; q=0.04_dp; r=0.08_dp; t=0.5_dp; sigma=0.25_dp
  call_res=european_option('call',spot,strike,q,r,t,sigma)
  out_res=barrier_option('downout','call',spot,strike,q,r,t,sigma,95.0_dp,0.0_dp,1400)
  in_res=barrier_option('downin','call',spot,strike,q,r,t,sigma,95.0_dp,0.0_dp,1400)
  call assert_close(in_res%value+out_res%value,call_res%value,5.0e-8_dp, &
    'Barrier in-out parity')
  call assert_true(out_res%value>=0.0_dp.and.in_res%value>=0.0_dp, &
    'Barrier values nonnegative')

  call european_implied_volatility('call',call_res%value,spot,strike,q,r,t,iv,status)
  call assert_true(status==0,'European implied-volatility status')
  call assert_close(iv,sigma,2.0e-8_dp,'European implied-volatility recovery')

  res=american_option('put',100.0_dp,105.0_dp,0.01_dp,0.04_dp,1.0_dp,0.25_dp,500)
  call american_implied_volatility('put',res%value,100.0_dp,105.0_dp, &
    0.01_dp,0.04_dp,1.0_dp,iv,status,500)
  call assert_true(status==0,'American implied-volatility status')
  call assert_close(iv,0.25_dp,2.0e-5_dp,'American implied-volatility recovery')

  call arithmetic_asian_mc('call',100.0_dp,100.0_dp,0.0_dp,0.03_dp, &
    1.0_dp,0.2_dp,12,20000,mc,se,1234)
  call assert_true(mc>0.0_dp.and.se<0.1_dp,'Arithmetic Asian Monte Carlo')

  call european_option_array('call',100.0_dp,[90.0_dp,100.0_dp,110.0_dp], &
    0.0_dp,0.03_dp,[0.5_dp,1.0_dp],0.2_dp,values)
  call assert_true(all(values(:,2)>=values(:,1)),'Option array maturity monotonicity')
  call assert_true(all(values(1,:)>=values(2,:)) .and. &
                   all(values(2,:)>=values(3,:)),'Option array strike monotonicity')
  write(*,'(a)') 'Option pricing and implied-volatility tests passed.'
end program test_options
