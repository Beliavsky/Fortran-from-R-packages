program test_blp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  integer, parameter :: jalt=2,t=2,k=1,h=40
  real(dp) :: r(1),sig(1,1),x(jalt*t,k),vdraw(k,h),share(jalt*t)
  real(dp), allocatable :: mu(:)
  real(dp), allocatable :: cp(:,:)
  real(dp) :: sigmar(1),a(1,1),th(1),cand(1,1),theta0(1),r0(1)
  type(blp_result) :: out
  integer :: i,info
  call rng_seed(808)
  r=0.0_dp; sig=r2sigma(r,1)
  if (abs(sig(1,1)-1.0_dp)>1.0e-14_dp) error stop "test_blp: r2sigma"
  x(:,1)=[0.0_dp,1.0_dp,0.4_dp,-0.5_dp]
  do i=1,h; vdraw(1,i)=0.5_dp*sin(real(i,dp)); end do
  share=[0.25_dp,0.35_dp,0.30_dp,0.20_dp]
  call share2mu(sig,x,vdraw,share,jalt,1.0e-8_dp,mu,cp,info,2000)
  if (info/=0 .or. .not.all(ieee_is_finite(mu))) error stop "test_blp: contraction"
  sigmar=1.0_dp; a=0.1_dp; th=0.0_dp; cand=0.05_dp; theta0=0.0_dp; r0=0.0_dp
  out=rbayes_blp(x,share,jalt,vdraw,sigmar,a,th,8,2,0.05_dp,cand,1.0e-7_dp,theta0,r0,1.0_dp, &
    nu0=4.0_dp,s0sq=1.0_dp)
  if (.not.all(ieee_is_finite(out%thetadraw)) .or. out%accept<0.0_dp .or. out%accept>1.0_dp) &
    error stop "test_blp: sampler"
  print *, "test_blp: PASS"
end program test_blp
