program test_diagnostics
   use compoissonreg
   implicit none
   integer,parameter :: n=20
   integer :: y(n),fails
   real(dp) :: x(n,1),s(n,1),res(n),lev(n),dev(n),boot(2,2)
   type(cmp_fit_t) :: fit,fit_free
   type(cmp_init_t) :: ini
   type(cmp_fixed_t) :: fix
   type(equitest_t) :: eq
   logical :: okboot(2)
   fails=0
   y=[0,1,2,1,3,2,4,1,0,2,1,2,3,1,0,4,2,1,3,2]
   x=1.0_dp;s=1.0_dp
   ini=default_init(1,1);fix=default_fixed(1,1);fix%gamma=.true.
   call fit_cmp_raw(y,x,s,fit,ini,fix)
   if(.not.(aic_cmp(fit)>0.0_dp.and.bic_cmp(fit)>0.0_dp))then
      print *,'FAIL AIC/BIC';fails=fails+1
   end if
   call residuals_cmp_raw(fit,res)
   if(any(.not.(abs(res)<huge(1.0_dp))))then;print *,'FAIL residuals';fails=fails+1;end if
   call bootstrap_cmp(fit,2,boot,okboot)
   if(any(.not.(abs(boot)<huge(1.0_dp))))then;print *,'FAIL bootstrap';fails=fails+1;end if

   ini=default_init(1,1);ini%beta=fit%beta;ini%gamma=0.0_dp
   call fit_cmp_raw(y,x,s,fit_free,ini)
   if(.not.fit_free%converged)then;print *,'FAIL free fit';fails=fails+1;end if
   eq=equitest_cmp(fit_free)
   if(eq%p_value<0.0_dp.or.eq%p_value>1.0_dp)then;print *,'FAIL equitest';fails=fails+1;end if
   call leverage_cmp(fit_free,lev)
   call deviance_cmp(fit_free,dev)
   if(any(lev<0.0_dp).or.any(.not.(abs(lev)<huge(1.0_dp))))then
      print *,'FAIL leverage';fails=fails+1
   end if
   if(any(dev<0.0_dp).or.any(.not.(abs(dev)<huge(1.0_dp))))then
      print *,'FAIL deviance';fails=fails+1
   end if
   if(fails==0)then;print *,'test_diagnostics: PASS';else;error stop 1;end if
end program test_diagnostics
