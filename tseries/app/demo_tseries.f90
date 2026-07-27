! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

program demo_tseries
   use tseries, only : dp, test_result, arma_result, garch_result, bds_result, &
      seed_random, quadratic_map, jarque_bera_test, adf_test, kpss_test, &
      arma_fit, garch_fit, bds_test, maximum_drawdown, sharpe_ratio, drawdown_result
   implicit none

   real(dp), allocatable :: x(:)
   type(test_result) :: jb, adf, kpss
   type(arma_result) :: ar
   type(garch_result) :: ga
   type(bds_result) :: bds
   type(drawdown_result) :: dd
   integer :: n

   n=400
   call seed_random(12345)
   x=quadratic_map(0.2_dp,3.7_dp,n)

   jb=jarque_bera_test(x)
   adf=adf_test(x,lags=2)
   kpss=kpss_test(x)
   ar=arma_fit(x,1,0,max_iterations=500)
   ga=garch_fit(x-sum(x)/real(n,dp),1,1,max_iterations=500)
   bds=bds_test(x,max_embedding=3,eps=[0.10_dp,0.20_dp])
   dd=maximum_drawdown(x)

   print '(a,i0)', 'observations: ',n
   print '(a,f12.6)', 'Jarque-Bera statistic: ',jb%statistic
   print '(a,f12.6)', 'ADF statistic: ',adf%statistic
   print '(a,f12.6)', 'KPSS statistic: ',kpss%statistic
   print '(a,*(f12.6,1x))', 'ARMA coefficients: ',ar%coefficients
   print '(a,*(f12.6,1x))', 'GARCH coefficients: ',ga%coefficients
   print '(a,f12.6)', 'BDS m=2, eps=0.1: ',bds%statistic(1,1)
   print '(a,f12.6)', 'maximum drawdown: ',dd%maximum
   print '(a,f12.6)', 'Sharpe ratio of changes: ',sharpe_ratio(x)
end program demo_tseries
