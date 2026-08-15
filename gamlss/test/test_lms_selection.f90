program test_lms_selection
   use gamlss
   implicit none
   integer,parameter::n=100
   real(dp)::x(n),y(n),prob(3)
   real(dp),allocatable::cent(:,:),mu(:),sig(:),nu(:)
   type(lms_result_t)::lm
   type(gamlss_family_comparison_t)::cmp
   type(gamlss_control_t)::ctl
   integer::i,status,best
   do i=1,n
      x(i)=real(i,dp)/real(n,dp)
      y(i)=10.0_dp+2.0_dp*x(i)+0.55_dp*sin(13.0_dp*real(i,dp))
   end do
   ctl%n_cyc=25;ctl%inner_cyc=25;ctl%c_crit=1.0e-5_dp
   call fit_lms(y,x,lm,family=GAMLSS_BCCG,df_mu=7,df_sigma=4,df_nu=4, &
      lambda_mu=4.0_dp,lambda_sigma=8.0_dp,lambda_nu=8.0_dp,control=ctl,status=status)
   call assert_true(status==0,'LMS fit')
   call predict_lms(lm,x,mu,sig,nu,status=status)
   call assert_true(status==0.and.all(mu>0.0_dp).and.all(sig>0.0_dp),'LMS prediction')
   prob=[0.1_dp,0.5_dp,0.9_dp]
   call lms_centiles(lm,x,prob,cent,status)
   call assert_true(status==0.and.all(cent(:,1)<cent(:,2)).and.all(cent(:,2)<cent(:,3)),'LMS centiles')

   call compare_families(y,[GAMLSS_NO,GAMLSS_TF],cmp,control=ctl)
   best=best_family(cmp)
   call assert_true(best==GAMLSS_NO.or.best==GAMLSS_TF,'family comparison')
   print '(a)','test_lms_selection: PASS'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then;print '(a)',trim(msg);error stop 1;end if
   end subroutine assert_true
end program test_lms_selection
