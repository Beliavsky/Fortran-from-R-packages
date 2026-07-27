! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program demo_rquantlib
  use rq_kinds, only: dp
  use rq_options
  use rq_curves
  use rq_bonds
  use rq_sabr
  use rq_hull_white
  implicit none
  type(option_result) :: euro, american, asian
  type(discount_curve_t) :: curve
  type(bond_result) :: bond
  type(sabr_result) :: sabr
  real(dp) :: strikes(5), vols(5), swaption
  integer :: status

  euro=european_option('call',100.0_dp,100.0_dp,0.01_dp,0.04_dp,1.0_dp,0.2_dp)
  american=american_option('put',100.0_dp,100.0_dp,0.01_dp,0.04_dp,1.0_dp,0.2_dp,500)
  asian=geometric_asian_option('call',100.0_dp,100.0_dp,0.01_dp,0.04_dp,1.0_dp,0.2_dp)
  write(*,'(a,f12.6)') 'European call:       ',euro%value
  write(*,'(a,f12.6)') 'American put:        ',american%value
  write(*,'(a,f12.6)') 'Geometric Asian call:',asian%value

  call make_flat_curve(0.04_dp,30.0_dp,0.25_dp,curve)
  call fixed_rate_bond_from_curve(100.0_dp,0.045_dp,5.0_dp,2,curve,bond)
  write(*,'(a,f12.6)') 'Fixed-rate bond NPV: ',bond%npv
  write(*,'(a,f12.6)') 'Bond yield:          ',bond%yield

  strikes=[0.02_dp,0.025_dp,0.03_dp,0.035_dp,0.04_dp]
  vols=sabr_lognormal_volatility(0.03_dp,strikes,5.0_dp,0.035_dp,0.5_dp,-0.2_dp,0.55_dp)
  call calibrate_sabr(0.03_dp,strikes,5.0_dp,vols,sabr,beta_fixed=0.5_dp)
  write(*,'(a,4f11.6)') 'SABR alpha,beta,rho,nu: ',sabr%alpha,sabr%beta,sabr%rho,sabr%nu

  call hull_white_swaption('payer',curve,0.1_dp,0.01_dp,2.0_dp, &
    [2.5_dp,3.0_dp,3.5_dp,4.0_dp],[0.5_dp,0.5_dp,0.5_dp,0.5_dp], &
    0.04_dp,1000000.0_dp,swaption,status)
  write(*,'(a,f14.4)') 'Hull-White payer swaption: ',swaption
end program demo_rquantlib
