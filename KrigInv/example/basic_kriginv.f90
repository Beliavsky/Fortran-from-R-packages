program basic_kriginv
  use kriginv
  implicit none
  type(krig_model) :: model
  type(optimization_result) :: opt
  type(optimizer_control) :: ctl
  real(dp) :: x(5,1),y(5),ell(1)
  logical :: ok
  x(:,1)=[0.0_dp,0.2_dp,0.5_dp,0.8_dp,1.0_dp]
  y=(x(:,1)-0.45_dp)**2-0.04_dp
  ell=0.3_dp
  call init_krig_model(model,x,y,ell,variance=0.2_dp,covariance='gauss',trend_order=0,ok=ok)
  if(.not.ok) error stop 'could not initialize kriging model'
  ctl%method='de'; ctl%pop_size=30; ctl%max_generations=25; ctl%seed=2026
  opt=max_infill_criterion([0.0_dp],[1.0_dp],'ranjan',[0.0_dp],model,1.0_dp,ctl)
  write(*,'(a,f10.6)') 'next point: ',opt%par(1,1)
  write(*,'(a,es12.4)') 'criterion:  ',opt%value
end program basic_kriginv
