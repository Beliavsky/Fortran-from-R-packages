! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program test_sabr_hull_white
  use rq_kinds, only: dp
  use rq_curves
  use rq_sabr
  use rq_hull_white
  use test_support
  implicit none
  type(sabr_result)::sabr_fit
  type(discount_curve_t)::curve
  type(hull_white_calibration_result)::hw_fit
  real(dp)::strikes(7),vols(7),forward,t,call,put,p0t,p0s
  real(dp)::payments(4),accruals(4),swaption,cap_prices(4)
  real(dp)::starts(4),ends(4),cap_strikes(4),notionals(4)
  integer::i,status

  forward=0.035_dp
  t=5.0_dp
  strikes=[0.015_dp,0.0225_dp,0.03_dp,0.035_dp,0.04_dp,0.05_dp,0.06_dp]
  vols=sabr_lognormal_volatility(forward,strikes,t,0.04_dp,0.5_dp,-0.25_dp,0.6_dp)
  call calibrate_sabr(forward,strikes,t,vols,sabr_fit,beta_fixed=0.5_dp)
  call assert_true(sabr_fit%status==0,'SABR calibration status')
  call assert_true(sabr_fit%rmse<2.0e-5_dp,'SABR synthetic calibration')
  call assert_true(sabr_fit%alpha>0.0_dp.and.sabr_fit%nu>0.0_dp, &
    'SABR parameter constraints')
  call assert_true(sabr_lognormal_volatility(forward,forward,t,0.04_dp, &
    0.5_dp,-0.25_dp,0.6_dp)>0.0_dp,'SABR ATM path')

  call make_flat_curve(0.03_dp,30.0_dp,0.1_dp,curve)
  p0t=curve%discount(5.0_dp)
  call assert_close(hull_white_discount_bond(curve,0.1_dp,0.01_dp,0.0_dp, &
    5.0_dp,0.03_dp),p0t,1.0e-12_dp,'Hull-White time-zero bond')
  p0s=curve%discount(2.0_dp)
  call=hull_white_bond_option('call',curve,0.1_dp,0.01_dp,2.0_dp,5.0_dp,0.8_dp)
  put=hull_white_bond_option('put',curve,0.1_dp,0.01_dp,2.0_dp,5.0_dp,0.8_dp)
  call assert_close(call-put,p0t-0.8_dp*p0s,1.0e-12_dp, &
    'Hull-White bond option parity')

  payments=[2.5_dp,3.0_dp,3.5_dp,4.0_dp]
  accruals=0.5_dp
  call hull_white_swaption('payer',curve,0.1_dp,0.01_dp,2.0_dp, &
    payments,accruals,0.035_dp,1000000.0_dp,swaption,status)
  call assert_true(status==0.and.swaption>=0.0_dp,'Hull-White swaption')

  starts=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  ends=starts+0.5_dp
  cap_strikes=0.035_dp
  notionals=1000000.0_dp
  do i=1,4
    cap_prices(i)=hull_white_caplet(curve,0.12_dp,0.015_dp,starts(i), &
      ends(i),cap_strikes(i),notionals(i))
  end do
  call calibrate_hull_white_caplets(curve,starts,ends,cap_strikes, &
    cap_prices,notionals,hw_fit)
  call assert_true(hw_fit%status==0,'Hull-White calibration status')
  call assert_true(hw_fit%rmse<1.0_dp,'Hull-White synthetic calibration')
  call assert_true(hw_fit%a>0.0_dp.and.hw_fit%sigma>0.0_dp, &
    'Hull-White positive parameters')
  write(*,'(a)') 'SABR and Hull-White tests passed.'
end program test_sabr_hull_white
