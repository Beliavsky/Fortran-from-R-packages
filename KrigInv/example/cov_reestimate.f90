program cov_reestimate
  use kriginv, only : dp, krig_model, km_control, fit_krig_model, update_krig_model
  implicit none
  integer, parameter :: n=8
  real(dp) :: x(n,1),y(n),xnew(1,1),ynew(1),old_range
  type(krig_model) :: model
  type(km_control) :: ctl
  integer :: i
  logical :: ok

  do i=1,n
    x(i,1)=real(i-1,dp)/real(n-1,dp)
    y(i)=sin(6.283185307179586_dp*x(i,1))+0.15_dp*x(i,1)
  end do

  ctl%pop_size=24
  ctl%multistart=3
  call fit_krig_model(model,x,y,covariance='matern5_2',trend_order=0,control=ctl,ok=ok)
  if(.not.ok) error stop 'initial kriging fit failed'

  old_range=model%lengthscale(1)
  xnew(1,1)=1.1_dp
  ynew(1)=sin(6.283185307179586_dp*xnew(1,1))+0.15_dp*xnew(1,1)

  ! Because this model estimated parameters, the default matches
  ! KrigInv/DiceKriging CovReEstimate=TRUE behavior.
  call update_krig_model(model,xnew,ynew,ok=ok)
  if(.not.ok) error stop 'kriging update failed'

  print '(a,f10.6)', 'range before update: ',old_range
  print '(a,f10.6)', 'range after update : ',model%lengthscale(1)
end program cov_reestimate
