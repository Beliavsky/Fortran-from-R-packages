program test_egi_reestimate
  use kriginv
  implicit none
  integer, parameter :: n=8
  real(dp) :: x(n,1),y(n),old_range
  type(krig_model) :: model
  type(egi_result) :: er_refit,er_fixed
  type(km_control) :: kc
  type(optimizer_control) :: oc
  integer :: i
  logical :: ok

  do i=1,n
    x(i,1)=real(i-1,dp)/real(n-1,dp)
    y(i)=fun(x(i,:))
  end do
  kc%pop_size=24; kc%multistart=3; kc%max_iter=250; kc%tol=1.0e-7_dp
  call seed_rng(314159)
  call fit_krig_model(model,x,y,covariance='matern5_2',trend_order=0,control=kc,ok=ok)
  if(.not.ok) error stop 'fitted model failed'
  old_range=model%lengthscale(1)

  oc%method='discrete'
  allocate(oc%optim_points(9,1))
  oc%optim_points(:,1)=[0.05_dp,0.15_dp,0.25_dp,0.35_dp,0.45_dp,0.55_dp,0.65_dp,0.85_dp,0.95_dp]

  call seed_rng(271828)
  er_refit=egi([0.0_dp],model,'ranjan',fun,1,[0.0_dp],[1.0_dp],opt_control=oc)
  if(.not.er_refit%ok) error stop 'default EGI failed'
  if(abs(er_refit%lastmodel%lengthscale(1)-old_range)<1.0e-5_dp) &
    error stop 'EGI default did not apply CovReEstimate'

  er_fixed=egi([0.0_dp],model,'ranjan',fun,1,[0.0_dp],[1.0_dp],opt_control=oc,cov_reestimate=.false.)
  if(.not.er_fixed%ok) error stop 'fixed EGI failed'
  if(abs(er_fixed%lastmodel%lengthscale(1)-old_range)>1.0e-12_dp) &
    error stop 'EGI CovReEstimate=FALSE changed range'

  print *, 'test_egi_reestimate: PASS'
contains
  real(dp) function fun(v) result(z)
    real(dp), intent(in) :: v(:)
    z=sin(6.283185307179586_dp*v(1))+0.15_dp*v(1)
  end function fun

  subroutine seed_rng(base)
    integer, intent(in) :: base
    integer :: ns,q
    integer, allocatable :: seed(:)
    call random_seed(size=ns); allocate(seed(ns))
    do q=1,ns
      seed(q)=base+32771*q
    end do
    call random_seed(put=seed)
  end subroutine seed_rng
end program test_egi_reestimate
