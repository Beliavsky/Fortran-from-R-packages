program test_batch
  use kriginv
  implicit none
  type(krig_model) :: model
  type(krig_prediction) :: pr
  type(integration_result) :: ir
  type(integration_control) :: ic
  type(optimization_result) :: opt
  type(optimizer_control) :: oc
  real(dp) :: x(6,1),y(6),ell(1),xn(1,1),s1,s2,s3,s4
  real(dp), allocatable :: pn(:),w(:)
  logical :: ok
  x(:,1)=[0.0_dp,0.15_dp,0.35_dp,0.60_dp,0.82_dp,1.0_dp]
  y=(x(:,1)-0.42_dp)**2-0.05_dp; ell=0.28_dp
  call init_krig_model(model,x,y,ell,variance=0.2_dp,nugget=0.0_dp,covariance='matern5_2',trend_order=1,ok=ok)
  call check(ok,'model')
  ic%n_points=30; ic%distrib='sobol'; ir=integration_design([0.0_dp],[1.0_dp],model,[0.0_dp],ic)
  pr=predict_nobias_km(model,ir%points,'UK',.false.)
  pn=excursion_probability(pr%mean,pr%sd,[0.0_dp])
  s1=sum(pn*(1.0_dp-pn))/real(size(pn),dp)
  xn(1,1)=0.48_dp
  s2=sur_optim_parallel(xn,ir%points,oldmean=pr%mean,oldsd=pr%sd,model=model,thresholds=[0.0_dp],current_sur=s1)
  call check(s2>=0.0_dp .and. s2<=s1+1.0e-8_dp,'SUR reduction')
  allocate(w(size(pr%sd))); w=1.0_dp
  s3=sum(pr%sd**2)/real(size(pr%sd),dp)
  s4=timse_optim_parallel(xn,ir%points,oldmean=pr%mean,oldsd=pr%sd,model=model,weight=w,current_timse=s3)
  call check(s4>=0.0_dp .and. s4<=s3+1.0e-8_dp,'IMSE reduction')
  oc%method='discrete'; allocate(oc%optim_points(5,1)); oc%optim_points(:,1)=[0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp]
  opt=max_infill_criterion([0.0_dp],[1.0_dp],'ranjan',[0.0_dp],model,1.0_dp,oc)
  call check(opt%ok .and. size(opt%par,1)==1,'discrete infill optimization')
  opt=max_sur_parallel([0.0_dp],[1.0_dp],1,ir,[0.0_dp],model,control=oc)
  call check(opt%ok .and. opt%value>=0.0_dp,'SUR optimization')
  opt=max_sur_parallel([0.0_dp],[1.0_dp],2,ir,[0.0_dp],model,control=oc)
  call check(opt%ok .and. size(opt%par,1)==2,'greedy discrete batch SUR')
  opt=max_timse_parallel([0.0_dp],[1.0_dp],2,ir,[0.0_dp],model,control=oc)
  call check(opt%ok .and. size(opt%par,1)==2,'greedy discrete batch TIMSE')
  opt=max_vorob_parallel([0.0_dp],[1.0_dp],2,ir,0.0_dp,model,control=oc)
  call check(opt%ok .and. size(opt%par,1)==2,'greedy discrete batch Vorob')
  ir%alpha=0.9_dp; ir%has_alpha=.true.
  opt=max_futurevol_parallel([0.0_dp],[1.0_dp],2,ir,0.0_dp,model,control=oc)
  call check(opt%ok .and. size(opt%par,1)==2,'greedy discrete future-volume')
  print '(a)', 'test_batch: PASS'
contains
  subroutine check(cond,msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if(.not.cond) then
      print '(a)', 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_batch
