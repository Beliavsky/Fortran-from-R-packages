program basic
   use gamlss
   implicit none
   integer,parameter::n=120
   real(dp)::x(n),y(n),xm(n,2),xs(n,1)
   type(gamlss_result_t)::fit
   type(gamlss_control_t)::ctl
   integer::i
   do i=1,n
      x(i)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=2.0_dp+1.4_dp*x(i)+0.35_dp*sin(11.0_dp*real(i,dp))
   end do
   xm(:,1)=1.0_dp;xm(:,2)=x;xs=1.0_dp
   ctl%n_cyc=30;ctl%c_crit=1.0e-6_dp
   call fit_gamlss_model(y,xm,GAMLSS_NO,fit,GAMLSS_METHOD_RS,x_sigma=xs,control=ctl)
   print '(a,l1)','Converged: ',fit%converged
   print '(a,2f10.4)','mu coefficients: ',fit%mu%coefficients
   print '(a,f10.4)','sigma: ',fit%sigma%fitted(1)
   print '(a,f12.4)','global deviance: ',fit%global_deviance
end program basic
