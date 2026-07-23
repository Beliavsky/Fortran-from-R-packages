! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

program test_arfima_backtests
   use rugarch
   implicit none
   type(arfima_spec) :: spec
   type(arfima_fit_result) :: fit
   type(var_test_result) :: vt
   type(directional_test_result) :: dt
   real(dp), allocatable :: x(:),fd(:),back(:),varf(:),actual(:),forecast(:)
   integer :: n

   call seed_rng(9876)
   spec=make_arfima_spec(1,0)
   spec%d=0.20_dp
   spec%ar(1)=0.25_dp
   allocate(x(900),fd(900),back(900))
   call simulate_arfima(spec,size(x),x)
   fit=fit_arfima(x,p=1)
   if (fit%status/=0) error stop 'ARFIMA fit failed'
   if (abs(fit%spec%d)>0.5_dp) error stop 'invalid fractional estimate'
   call fractional_difference(x,spec%d,fd)
   call fractional_integrate(fd,spec%d,back)
   if (sum(abs(back(200:)-x(200:)))/real(size(x)-199,dp)>0.20_dp) error stop 'fractional round trip failed'

   n=200
   allocate(actual(n),forecast(n),varf(n))
   actual=x(1:n)
   forecast=0.2_dp*x(2:n+1)
   varf=-1.65_dp
   vt=var_test(0.05_dp,actual,varf)
   if (vt%actual_exceedances<0) error stop 'VaR test failed'
   dt=directional_accuracy_test(forecast,actual)
   if (dt%directional_accuracy<0.0_dp .or. dt%directional_accuracy>1.0_dp) error stop 'DAC failed'
   print '(a)', 'ARFIMA and backtest tests passed'
end program test_arfima_backtests
