program multi_step_forecast
  use rumidas
  implicit none
  real(dp)::param(7)
  real(dp),allocatable::forecast(:)
  type(garch_midas_spec)::spec
  integer::i,status

  param=[0.06_dp,0.85_dp,0.04_dp,0.08_dp,-8.0_dp,0.2_dp,2.5_dp]
  spec=garch_midas_spec(RUMIDAS_GMX,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,3,0,.true.)
  call multi_step_ahead_pred(param,spec,1.25_dp,1.2e-4_dp,10,forecast,status,x_last=0.4_dp,x_ar1=0.7_dp)
  if(status/=0) error stop 'forecast failed'
  do i=1,size(forecast)
    print '(i3,es16.7)',i,forecast(i)
  end do
end program multi_step_forecast
