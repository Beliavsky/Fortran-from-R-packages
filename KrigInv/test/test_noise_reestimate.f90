program test_noise_reestimate
  use kriginv, only : dp, krig_model, krig_prediction, km_control, fit_krig_model, &
                      update_krig_model, predict_nobias_km
  implicit none
  integer, parameter :: n=9
  real(dp) :: x(n,1),y(n),xn(1,1),yn(1),nv(1),old_range
  type(krig_model) :: model
  type(krig_prediction) :: pred
  type(km_control) :: ctl
  integer :: i
  logical :: ok

  do i=1,n
    x(i,1)=real(i-1,dp)/real(n-1,dp)
    y(i)=cos(5.0_dp*x(i,1))+0.2_dp*x(i,1)
  end do
  ctl%pop_size=24; ctl%multistart=3; ctl%max_iter=250; ctl%tol=1.0e-7_dp
  call seed_rng(12345)
  call fit_krig_model(model,x,y,covariance='matern3_2',trend_order=1,control=ctl,ok=ok)
  if(.not.ok) error stop 'initial fit failed'
  old_range=model%lengthscale(1)

  xn(1,1)=0.43_dp; yn(1)=cos(5.0_dp*xn(1,1))+0.2_dp*xn(1,1); nv(1)=0.0025_dp
  call seed_rng(54321)
  call update_krig_model(model,xn,yn,nv,ok,cov_reestimate=.true.)
  if(.not.ok) error stop 'noisy covariance re-estimation failed'
  if(.not.model%dice_model%noise_flag) error stop 'new noise did not activate noise model'
  if(size(model%noise)/=n+1) error stop 'noise vector length mismatch'
  if(abs(model%noise(n+1)-nv(1))>1.0e-15_dp) error stop 'new noise variance mismatch'
  if(abs(model%lengthscale(1)-old_range)<1.0e-6_dp) error stop 'noisy update did not refit range'

  pred=predict_nobias_km(model,reshape([0.43_dp,0.8_dp],[2,1]),'UK',.true.)
  if(.not.pred%ok .or. any(pred%sd<0.0_dp)) error stop 'noisy prediction failed'
  if(any(.not.(pred%mean<huge(1.0_dp)))) error stop 'non-finite noisy prediction'

  print *, 'test_noise_reestimate: PASS'
contains
  subroutine seed_rng(base)
    integer, intent(in) :: base
    integer :: ns,q
    integer, allocatable :: seed(:)
    call random_seed(size=ns); allocate(seed(ns))
    do q=1,ns
      seed(q)=base+65537*q
    end do
    call random_seed(put=seed)
  end subroutine seed_rng
end program test_noise_reestimate
