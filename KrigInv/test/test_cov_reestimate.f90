program test_cov_reestimate
  use kriginv, only : dp, krig_model, km_control, fit_krig_model, update_krig_model
  implicit none
  integer, parameter :: n=16
  real(dp) :: x(n,2),y(n),xnew(1,2),ynew(1),old_range(2),old_var
  type(krig_model) :: model,refit,fixed
  type(km_control) :: ctl
  integer :: i,j,k
  logical :: ok

  k=0
  do j=0,3
    do i=0,3
      k=k+1
      x(k,1)=real(i,dp)/3.0_dp
      x(k,2)=real(j,dp)/3.0_dp
      y(k)=branin_local(x(k,:))
    end do
  end do
  ctl%pop_size=80; ctl%multistart=8; ctl%max_iter=500; ctl%tol=1.0e-8_dp
  call seed_rng(19081)
  call fit_krig_model(model,x,y,covariance='matern5_2',trend_order=0,control=ctl,ok=ok)
  if(.not.ok) error stop 'fit_krig_model failed'
  if(.not.model%uses_dicekriging) error stop 'DiceKriging backend not enabled'
  if(.not.model%param_estim .or. .not.model%cov_reestimate_default) error stop 'param.estim default mismatch'
  if(maxval(abs(model%lengthscale-[0.8254355_dp,2.0_dp]))>2.0e-3_dp) error stop 'Branin range mismatch'
  if(abs(model%variance-145556.6_dp)/145556.6_dp>5.0e-3_dp) error stop 'Branin variance mismatch'
  if(abs(model%beta(1)-306.5783_dp)>0.3_dp) error stop 'Branin trend mismatch'

  old_range=model%lengthscale; old_var=model%variance
  xnew(1,:)=[0.45_dp,0.55_dp]; ynew(1)=branin_local(xnew(1,:))

  refit=model
  call seed_rng(271828)
  call update_krig_model(refit,xnew,ynew,ok=ok)
  if(.not.ok) error stop 'default re-estimating update failed'
  if(refit%n/=n+1) error stop 're-estimating update did not append point'
  if(maxval(abs(refit%lengthscale-old_range))<1.0e-5_dp .and. &
     abs(refit%variance-old_var)<1.0e-3_dp) error stop 'default CovReEstimate did not refit parameters'

  fixed=model
  call update_krig_model(fixed,xnew,ynew,ok=ok,cov_reestimate=.false.)
  if(.not.ok) error stop 'fixed-covariance update failed'
  if(maxval(abs(fixed%lengthscale-old_range))>1.0e-12_dp) error stop 'CovReEstimate=FALSE changed range'
  if(abs(fixed%variance-old_var)>1.0e-8_dp) error stop 'CovReEstimate=FALSE changed variance'

  print *, 'test_cov_reestimate: PASS'
contains
  real(dp) function branin_local(xx) result(v)
    real(dp), intent(in) :: xx(:)
    real(dp) :: x1,x2,a,b,c,r,s,t,pi
    pi=acos(-1.0_dp); x1=15.0_dp*xx(1)-5.0_dp; x2=15.0_dp*xx(2)
    a=1.0_dp; b=5.0_dp/(4.0_dp*pi*pi); c=5.0_dp/pi
    r=6.0_dp; s=10.0_dp; t=1.0_dp/(8.0_dp*pi)
    v=a*(x2-b*x1*x1+c*x1-r)**2+s*(1.0_dp-t)*cos(x1)+s
  end function branin_local

  subroutine seed_rng(base)
    integer, intent(in) :: base
    integer :: ns,q
    integer, allocatable :: seed(:)
    call random_seed(size=ns); allocate(seed(ns))
    do q=1,ns
      seed(q)=base+104729*q
    end do
    call random_seed(put=seed)
  end subroutine seed_rng
end program test_cov_reestimate
