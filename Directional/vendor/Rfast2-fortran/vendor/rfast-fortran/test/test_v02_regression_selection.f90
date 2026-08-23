program test_v02_regression_selection
   use rfast
   implicit none
   integer,parameter::n=60
   real(dp)::x(n,3),y(n),yg(n),yp(n),yi(n),ymat(n,2)
   integer::yc(n),i,fail
   type(regression_result)::r
   type(selection_result)::s
   type(multinomial_result)::mr
   type(multivariate_regression_result)::sm
   fail=0
   do i=1,n
      x(i,1)=(real(i,dp)-30.5_dp)/15.0_dp;x(i,2)=sin(1.7_dp*real(i,dp));x(i,3)=cos(0.43_dp*real(i,dp))
      y(i)=1.5_dp+2.2_dp*x(i,1)+0.08_dp*sin(real(i,dp))
      yg(i)=exp(0.3_dp+0.45_dp*x(i,1))*(1.0_dp+0.04_dp*sin(real(i,dp)))
      yi(i)=exp(0.2_dp+0.25_dp*x(i,1))*(1.0_dp+0.03_dp*cos(real(i,dp)))
      yp(i)=real(mod(i,5),dp)+exp(0.2_dp+0.25_dp*x(i,1))
      yc(i)=mod(i-1,3)+1
      ymat(i,1)=1.0_dp+1.2_dp*x(i,1)+0.05_dp*sin(real(i,dp));ymat(i,2)=-0.5_dp+0.7_dp*x(i,1)+0.04_dp*cos(real(i,dp))
   end do
   s=ompr(y,x,OMP_BIC,2.0_dp,.true.);if(s%status/=0.or.s%selected(1)/=1)fail=fail+1
   s=bic_corfsreg(y,x,2.0_dp);if(s%status/=0.or.s%selected(1)/=1)fail=fail+1
   r=gamma_regression(yg,x(:,1:1));if(r%status/=0.or.abs(r%beta(2)-0.45_dp)>0.08_dp)fail=fail+1
   r=invgauss_regression(yi,x(:,1:1));if(r%status/=0.or.any(.not.(r%fitted>0.0_dp)))fail=fail+1
   r=quasipoisson_regression(yp,x(:,1:1));if(r%status/=0.or.r%dispersion<=0.0_dp)fail=fail+1
   r=proportion_regression(merge(1.0_dp,0.0_dp,x(:,1)>0.0_dp),x(:,1:1),.true.);if(r%status/=0)fail=fail+1
   mr=multinomial_regression(yc,x(:,1:2));if(mr%status/=0.or..not.allocated(mr%beta))fail=fail+1
   sm=spatial_median_regression(ymat,x(:,1:1));if(sm%status/=0.or.abs(sm%beta(2,1)-1.2_dp)>0.1_dp)fail=fail+1
   s=omp_glm(merge(1.0_dp,0.0_dp,x(:,1)>0.0_dp),x,OMP_LOGISTIC,2.0_dp,.true.)
   if(s%status/=0.or.s%selected(1)/=1)fail=fail+1
   if(fail==0)then
      print '(a)','test_v02_regression_selection: PASS'
   else
      print '(a,i0)','test_v02_regression_selection: FAIL ',fail
      error stop 1
   end if
end program test_v02_regression_selection
