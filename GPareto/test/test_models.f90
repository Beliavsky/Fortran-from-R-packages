program test_models
  use gpareto, only : dp, gp_model_set, fit_gp_model, trend_const, predict_gps
  use gpareto, only : nondominated_points, crit_ehi, crit_qehi, crit_sur, update_gp
  implicit none
  type(gp_model_set) :: ms
  real(dp) :: x(7,1), y1(7), y2(7), q(2,1), integ(9,1), ref(2), qv
  real(dp), allocatable :: mu(:,:),sd(:,:),obs(:,:),front(:,:),v(:)
  integer :: i
  do i=1,7;x(i,1)=real(i-1,dp)/6.0_dp;y1(i)=x(i,1)**2;y2(i)=(1.0_dp-x(i,1))**2;end do
  allocate(ms%model(2))
  call fit_gp_model(ms%model(1),x,y1,covtype='gauss',trend_kind=trend_const,nugget=1.0e-8_dp)
  call fit_gp_model(ms%model(2),x,y2,covtype='gauss',trend_kind=trend_const,nugget=1.0e-8_dp)
  q(:,1)=[0.25_dp,0.75_dp];call predict_gps(ms,q,mu,sd)
  if(any(.not.(abs(mu)<huge(1.0_dp))).or.any(sd<0.0_dp)) error stop 'prediction invalid'
  allocate(obs(7,2));obs(:,1)=y1;obs(:,2)=y2;call nondominated_points(obs,front);ref=[1.2_dp,1.2_dp]
  call crit_ehi(q,ms,front,ref,v,nsamp=100)
  if(any(v<0.0_dp)) error stop 'EHI invalid'
  call crit_qehi(q,ms,front,ref,qv,nsamp=20)
  if(qv<0.0_dp) error stop 'qEHI invalid'
  do i=1,9;integ(i,1)=real(i,dp)/10.0_dp;end do
  call crit_sur(q,ms,front,integ,v,nsamp=8)
  if(any(.not.(abs(v)<huge(1.0_dp)))) error stop 'SUR invalid'
  call update_gp(ms%model(1),reshape([0.37_dp],[1,1]),[0.37_dp**2],cov_reestimate=.false.)
  if(ms%model(1)%km%n/=8) error stop 'update failed'
  print *, 'test_models PASS'
end program test_models
