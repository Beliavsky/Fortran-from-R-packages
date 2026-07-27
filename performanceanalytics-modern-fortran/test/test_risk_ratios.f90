! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program test_risk_ratios
  use kinds_mod, only: dp
  use risk_mod
  use performance_ratios_mod
  use test_support_mod
  implicit none
  real(dp) :: r(5), rf(5), b(5), vars(5), comp(2), assets(5,2), w(2), total
  type(var_backtest_result) :: bt
  integer :: i

  r=[-0.10_dp,-0.05_dp,0.00_dp,0.02_dp,0.03_dp]
  rf=0.0_dp
  call assert_close(historical_var(r,0.8_dp),0.06_dp,1.0e-13_dp,'historical var')
  call assert_close(historical_es(r,0.8_dp),0.10_dp,1.0e-13_dp,'historical es')
  call assert_true(gaussian_es(r,0.95_dp)>gaussian_var(r,0.95_dp),'gaussian es greater var')
  call assert_true(modified_es(r,0.95_dp)>0.0_dp,'modified es finite')
  call assert_close(sample_var(r,0.2_dp,.true.),-0.06_dp,1.0e-13_dp,'sample var')
  call assert_true(rachev_tail_ratio(r,0.2_dp,0.2_dp)>0.0_dp,'rachev ratio')

  vars=0.04_dp
  call var_backtest(r,vars,0.8_dp,bt)
  call assert_true(bt%violations==2,'var violations')
  call assert_true(bt%kupiec_pvalue>=0.0_dp .and. bt%kupiec_pvalue<=1.0_dp,'kupiec p')

  assets(:,1)=r
  assets(:,2)=[0.01_dp,-0.01_dp,0.02_dp,0.00_dp,0.01_dp]
  w=[0.7_dp,0.3_dp]
  call component_var(assets,w,0.95_dp,'gaussian',comp,total)
  call assert_true(total>0.0_dp,'portfolio var total')
  call assert_true(all(abs(comp)<100.0_dp),'component var finite')

  call assert_close(sharpe_ratio([0.01_dp,0.02_dp,0.03_dp], [0.0_dp,0.0_dp,0.0_dp]), &
                    0.02_dp/0.01_dp,1.0e-12_dp,'sharpe ratio')
  call assert_true(sortino_ratio([0.01_dp,-0.02_dp,0.03_dp],0.0_dp)>0.0_dp,'sortino')
  call assert_true(omega_ratio([0.02_dp,-0.01_dp,0.03_dp,-0.02_dp],0.0_dp)>1.0_dp,'omega')
  call assert_true(kappa_ratio([0.02_dp,-0.01_dp,0.03_dp,-0.02_dp],0.0_dp,3)>0.0_dp,'kappa')
  call assert_true(bernardo_ledoit_ratio([0.02_dp,-0.01_dp,0.03_dp,-0.02_dp])>1.0_dp,'bernardo ledoit')
  call assert_true(calmar_ratio([0.02_dp,-0.01_dp,0.03_dp,-0.02_dp],12.0_dp)>0.0_dp,'calmar')
  call assert_true(pain_ratio([0.02_dp,-0.01_dp,0.03_dp,-0.02_dp],12.0_dp,0.0_dp)>0.0_dp,'pain ratio')
  b=[-0.08_dp,-0.03_dp,0.01_dp,0.01_dp,0.02_dp]
  call assert_true(tracking_error(r,b,12.0_dp)>0.0_dp,'tracking error')
  call assert_true(probabilistic_sharpe_ratio(r,rf,0.0_dp)>=0.0_dp,'probabilistic sharpe')
  call assert_true(minimum_track_record_length(r,rf,-1.0_dp,0.95_dp)>1.0_dp,'minimum track record')
  call assert_true(hurst_index([(0.001_dp*real(i,dp),i=1,128)])>=0.0_dp,'hurst')
  write(*,'(a)')'Risk and performance-ratio tests passed.'
end program test_risk_ratios
