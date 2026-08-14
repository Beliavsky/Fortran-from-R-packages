program test_core
  use kriginv
  use kriginv_math, only : bvn_cdf
  implicit none
  type(krig_model) :: model
  type(krig_prediction) :: pr
  type(integration_result) :: ir
  type(integration_control) :: ic
  real(dp) :: x(5,1),y(5),ell(1),q
  real(dp), allocatable :: c(:),ep(:)
  logical :: ok
  x(:,1)=[0.0_dp,0.25_dp,0.5_dp,0.75_dp,1.0_dp]
  y=sin(6.283185307179586_dp*x(:,1)); ell=0.35_dp
  call init_krig_model(model,x,y,ell,variance=1.0_dp,nugget=0.0_dp,covariance='gauss',trend_order=0,ok=ok)
  call check(ok,'model initialization')
  pr=predict_nobias_km(model,x,'UK',.true.)
  call check(maxval(abs(pr%mean-y))<1.0e-8_dp,'interpolation mean')
  call check(maxval(pr%sd)<1.0e-6_dp,'interpolation sd')
  call check(maxval(abs(pr%covariance-transpose(pr%covariance)))<1.0e-12_dp,'prediction covariance symmetry')
  ep=excursion_probability([0.0_dp,1.0_dp], [1.0_dp,1.0_dp], [0.0_dp])
  call check(abs(ep(1)-0.5_dp)<1.0e-14_dp,'excursion p=0.5')
  q=bvn_cdf(0.0_dp,0.0_dp,0.5_dp)
  call check(abs(q-(0.25_dp+asin(0.5_dp)/(2.0_dp*acos(-1.0_dp))))<2.0e-10_dp,'bivariate normal cdf')
  ic%n_points=20; ic%distrib='sobol'
  ir=integration_design([0.0_dp],[1.0_dp],model,[0.0_dp],ic)
  call check(ir%ok .and. size(ir%points,1)==20,'Sobol integration design')
  c=ranjan_optim(ir%points,model,0.0_dp,1.0_dp)
  call check(all(c>=-1.0e-14_dp),'Ranjan criterion finite/nonnegative')
  c=bichon_optim(ir%points,model,0.0_dp,1.0_dp)
  call check(all(c>=-1.0e-14_dp),'Bichon criterion finite/nonnegative')
  c=tmse_optim(ir%points,model,[0.0_dp],0.0_dp)
  call check(all(c>=-1.0e-14_dp),'TMSE criterion finite/nonnegative')
  c=tsee_optim(ir%points,model,0.0_dp)
  call check(all(c>=-1.0e-14_dp),'TSEE criterion finite/nonnegative')
  print '(a)', 'test_core: PASS'
contains
  subroutine check(cond,msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if(.not.cond) then
      print '(a)', 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_core
