program test_core
   use gamlss
   implicit none
   integer,parameter::n=160
   real(dp)::xv(n),y(n),xmu(n,2),xs(n,1),mse,m,sd,sk,ku
   real(dp),allocatable::rq(:),basis(:,:),pen(:,:),basis2(:,:)
   type(gamlss_result_t)::rs,cg,sm
   type(gamlss_control_t)::ctl
   type(p_spline_spec_t)::spec
   integer::i,status

   do i=1,n
      xv(i)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.25_dp+1.8_dp*xv(i)+0.42_dp*sin(17.0_dp*real(i,dp))
   end do
   xmu(:,1)=1.0_dp;xmu(:,2)=xv;xs=1.0_dp
   ctl%n_cyc=35;ctl%inner_cyc=30;ctl%c_crit=1.0e-6_dp;ctl%inner_crit=1.0e-7_dp
   call fit_gamlss_model(y,xmu,GAMLSS_NO,rs,GAMLSS_METHOD_RS,x_sigma=xs,control=ctl)
   call assert_true(rs%status==0,'RS status')
   call assert_true(abs(rs%mu%coefficients(1)-1.25_dp)<0.08_dp,'RS intercept')
   call assert_true(abs(rs%mu%coefficients(2)-1.8_dp)<0.08_dp,'RS slope')
   call assert_true(rs%sigma%fitted(1)>0.2_dp.and.rs%sigma%fitted(1)<0.7_dp,'RS sigma')

   ctl%n_cyc=25
   call fit_gamlss_model(y,xmu,GAMLSS_NO,cg,GAMLSS_METHOD_CG,x_sigma=xs,control=ctl)
   call assert_true(cg%status==0,'CG status')
   call assert_true(abs(cg%mu%coefficients(2)-1.8_dp)<0.1_dp,'CG slope')

   call randomized_quantile_residuals(y,rs,rq,status,randomize=.false.)
   call assert_true(status==0.and.all(abs(rq)<10.0_dp),'quantile residuals')
   call residual_moments(rq,m,sd,sk,ku)
   call assert_true(sd>0.1_dp,'residual moments')

   call fit_p_spline_basis(xv,spec,basis,df=9,status=status)
   call assert_true(status==0.and.size(basis,1)==n,'P-spline basis')
   call assert_true(maxval(abs(spec%penalty-transpose(spec%penalty)))<1.0e-12_dp,'penalty symmetry')
   call predict_p_spline_basis(xv,spec,basis2,status)
   call assert_true(status==0.and.maxval(abs(basis-basis2))<1.0e-10_dp,'persistent spline basis')
   pen=spec%penalty
   y=2.0_dp+0.7_dp*sin(3.141592653589793_dp*xv)+0.08_dp*cos(19.0_dp*real([(i,i=1,n)],dp))
   call fit_gamlss_model(y,basis,GAMLSS_NO,sm,GAMLSS_METHOD_RS,x_sigma=xs, &
      penalty_mu=pen,lambda_mu=2.0_dp,control=ctl)
   call assert_true(sm%status==0,'smooth fit status')
   mse=sum((sm%mu%fitted-(2.0_dp+0.7_dp*sin(3.141592653589793_dp*xv)))**2)/real(n,dp)
   call assert_true(mse<0.03_dp,'smooth mean recovery')

   print '(a)','test_core: PASS'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then;print '(a)',trim(msg);error stop 1;end if
   end subroutine assert_true
end program test_core
