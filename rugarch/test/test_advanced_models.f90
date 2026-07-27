! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

program test_advanced_models
   use rugarch
   implicit none
   integer,parameter :: n=140
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   real(dp) :: y(n),sigma(n),residuals(n),realized(n),measurement_residuals(n)
   real(dp) :: weights(60)
   logical :: valid

   call seed_rng(77123)

   spec=make_garch_spec(2,1,model_figarch,dist_norm)
   spec%omega=1.0e-5_dp
   spec%alpha=[0.10_dp,0.04_dp]
   spec%beta=[0.35_dp]
   spec%frac_d=0.30_dp
   spec%figarch_truncation=60
   call figarch_weights(spec%frac_d,spec%alpha,spec%beta,weights)
   if(any(weights<0.0_dp)) error stop 'invalid generalized FIGARCH weights'
   call simulate_garch(spec,n,y,sigma,residuals,burn_in=100)
   call garch_filter(y,spec,residuals,sigma,valid)
   if(.not.valid) error stop 'FIGARCH filter failed'
   fit=fit_figarch(y,p=2,q=1,truncation=60,max_iterations=100)
   call check_fit(fit,'FIGARCH')

   spec=make_garch_spec(1,1,model_csgarch,dist_norm)
   spec%omega=1.0e-5_dp
   spec%alpha=0.05_dp
   spec%beta=0.70_dp
   spec%rho=0.95_dp
   spec%phi=0.03_dp
   call simulate_garch(spec,n,y,sigma,residuals,burn_in=100)
   fit=fit_csgarch11(y,max_iterations=100)
   call check_fit(fit,'component GARCH')

   spec=make_garch_spec(1,1,model_realgarch,dist_norm)
   spec%omega=-0.20_dp
   spec%alpha=0.15_dp
   spec%beta=0.70_dp
   spec%xi=-0.10_dp
   spec%phi=0.90_dp
   spec%tau1=-0.10_dp
   spec%tau2=0.05_dp
   spec%measurement_sd=0.15_dp
   call simulate_realgarch(spec,n,y,realized,sigma,residuals,burn_in=100)
   call realgarch_filter(y,realized,spec,residuals,sigma,measurement_residuals,valid)
   if(.not.valid) error stop 'realGARCH filter failed'
   fit=fit_realgarch11(y,realized,max_iterations=120)
   call check_fit(fit,'realGARCH')
   if(.not.allocated(fit%measurement_residuals)) error stop 'realGARCH measurement residuals missing'

   spec=make_garch_spec(1,1,model_fgarch,dist_norm)
   call configure_fgarch_submodel(spec,fgarch_allgarch)
   spec%omega=1.0e-5_dp
   spec%alpha=0.08_dp
   spec%beta=0.82_dp
   spec%eta1=0.10_dp
   spec%eta2=0.15_dp
   spec%fgarch_lambda=1.50_dp
   call simulate_garch(spec,n,y,sigma,residuals,burn_in=100)
   call garch_filter(y,spec,residuals,sigma,valid)
   if(.not.valid) error stop 'fGARCH filter failed'
   fit=fit_fgarch11(y,max_iterations=120,submodel=fgarch_allgarch)
   call check_fit(fit,'fGARCH')
   if(fit%spec%fgarch_submodel/=fgarch_allgarch) error stop 'fGARCH submodel not retained'

   call check_submodels
   print '(a)', 'advanced model tests passed'

contains
   subroutine check_fit(result,label)
      type(garch_fit_result),intent(in)::result
      character(len=*),intent(in)::label
      if(.not.allocated(result%sigma)) error stop label//' fit output missing'
      if(any(result%sigma<=0.0_dp)) error stop label//' fit has nonpositive sigma'
      if(.not.(result%log_likelihood>-huge(1.0_dp)/10.0_dp)) error stop label//' fit likelihood invalid'
   end subroutine check_fit

   subroutine check_submodels
      type(garch_spec)::s
      type(garch_fit_result)::subfit
      integer::submodel
      do submodel=fgarch_garch,fgarch_gjrgarch
         s=make_garch_spec(1,1,model_fgarch,dist_norm)
         call configure_fgarch_submodel(s,submodel)
         if(s%fgarch_lambda<=0.0_dp) error stop 'invalid fGARCH submodel power'
         if(s%delta+s%fgarch_fk*s%fgarch_lambda<=0.0_dp) error stop 'invalid fGARCH shock power'
         subfit=fit_fgarch11(y,fit_mean=.false.,max_iterations=5,submodel=submodel)
         if(.not.allocated(subfit%sigma)) error stop 'fGARCH submodel fitter failed'
         if(subfit%spec%fgarch_submodel/=submodel) error stop 'fGARCH submodel identifier lost'
      end do
   end subroutine check_submodels
end program test_advanced_models
