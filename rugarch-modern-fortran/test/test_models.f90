! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

program test_models
   use rugarch
   implicit none
   integer, parameter :: n=700
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   real(dp), allocatable :: y(:),s(:),e(:),ef(:),sf(:),fc(:)
   logical :: ok

   call seed_rng(12345)
   spec=make_garch_spec(1,1,model_gjrgarch,dist_norm)
   spec%omega=1.0e-5_dp
   spec%alpha=0.07_dp
   spec%gamma=0.08_dp
   spec%beta=0.85_dp
   allocate(y(n),s(n),e(n),ef(n),sf(n),fc(3))
   call simulate_garch(spec,n,y,s,e)
   call garch_filter(y,spec,ef,sf,ok)
   if (.not.ok) error stop 'GJR filter failed'
   if (any(sf<=0.0_dp)) error stop 'nonpositive sigma'
   if (.not.(garch_log_likelihood(y,spec)>-huge(1.0_dp)/10.0_dp)) error stop 'invalid likelihood'

   fit=fit_garch11(y,dist_norm,max_iterations=700)
   if (.not.allocated(fit%sigma)) error stop 'fit output missing'
   if (fit%spec%alpha(1)<0.0_dp .or. fit%spec%beta(1)<0.0_dp) error stop 'invalid fit'
   call forecast_volatility(fit%spec,fit%residuals,fit%sigma,3,fc)
   if (any(fc<=0.0_dp)) error stop 'invalid forecast'
   print '(a)', 'model tests passed'
end program test_models
