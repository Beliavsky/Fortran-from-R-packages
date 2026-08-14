program basic_gpareto
  use gpareto, only : dp, gp_model_set, fit_gp_model, trend_const
  use gpareto, only : nondominated_points, crit_ehi
  implicit none
  type(gp_model_set) :: models
  real(dp) :: x(7,1), y1(7), y2(7), q(1,1), ref(2)
  real(dp), allocatable :: obs(:,:), front(:,:), ehi(:)
  integer :: i

  do i=1,7
    x(i,1)=real(i-1,dp)/6.0_dp
    y1(i)=x(i,1)**2
    y2(i)=(1.0_dp-x(i,1))**2
  end do

  allocate(models%model(2))
  call fit_gp_model(models%model(1),x,y1,covtype='gauss', &
    trend_kind=trend_const,nugget=1.0e-8_dp)
  call fit_gp_model(models%model(2),x,y2,covtype='gauss', &
    trend_kind=trend_const,nugget=1.0e-8_dp)

  allocate(obs(7,2))
  obs(:,1)=y1
  obs(:,2)=y2
  call nondominated_points(obs,front)
  q(1,1)=0.42_dp
  ref=[1.2_dp,1.2_dp]
  call crit_ehi(q,models,front,ref,ehi)
  print '(a,f12.8)', 'EHI at x=0.42: ',ehi(1)
end program basic_gpareto
