program test_egi
  use kriginv
  implicit none
  type(krig_model) :: model
  type(egi_result) :: er
  type(optimizer_control) :: oc
  type(integration_control) :: ic
  real(dp) :: x(4,1),y(4),ell(1)
  logical :: ok
  x(:,1)=[0.0_dp,0.30_dp,0.70_dp,1.0_dp]
  y=[f(x(1,:)),f(x(2,:)),f(x(3,:)),f(x(4,:))]
  ell=0.30_dp
  call init_krig_model(model,x,y,ell,variance=0.2_dp,covariance='gauss',trend_order=0,ok=ok)
  call check(ok,'initial model')
  oc%method='discrete'; allocate(oc%optim_points(7,1)); oc%optim_points(:,1)=[0.1_dp,0.2_dp,0.4_dp,0.5_dp,0.6_dp,0.8_dp,0.9_dp]
  er=egi([0.0_dp],model,'ranjan',f,2,[0.0_dp],[1.0_dp],opt_control=oc)
  call check(er%ok .and. er%nsteps==2 .and. er%lastmodel%n==6,'sequential EGI')
  ic%n_points=20; ic%distrib='sobol'
  er=egi_parallel([0.0_dp],model,'sur',f,1,1,[0.0_dp],[1.0_dp],opt_control=oc,int_control=ic)
  call check(er%ok .and. er%nsteps==1 .and. er%lastmodel%n==5,'SUR EGI')
  print '(a)', 'test_egi: PASS'
contains
  real(dp) function f(v) result(z)
    real(dp), intent(in) :: v(:)
    z=(v(1)-0.45_dp)**2-0.04_dp
  end function f
  subroutine check(cond,msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if(.not.cond) then
      print '(a)', 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_egi
