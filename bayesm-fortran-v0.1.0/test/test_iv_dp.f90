program test_iv_dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  integer, parameter :: n=20
  real(dp) :: y(n),x(n),z(n,2),w(n,1),mbg(2),abg(2,2),md(2),ad(2,2),v2(2,2)
  real(dp) :: yd(12,1),m0(1),s0(1,1)
  type(iv_result) :: ivr
  type(dp_mixture_result) :: dpr
  integer :: i
  call rng_seed(2026)
  do i=1,n
    z(i,1)=1.0_dp; z(i,2)=real(i-n/2,dp)/10.0_dp; w(i,1)=1.0_dp
    x(i)=0.5_dp+0.8_dp*z(i,2)+0.1_dp*sin(real(i,dp)); y(i)=1.0_dp+1.5_dp*x(i)+0.05_dp*cos(real(i,dp))
  end do
  mbg=0.0_dp; abg=0.0_dp; abg(1,1)=0.01_dp; abg(2,2)=0.01_dp; md=0.0_dp; ad=0.0_dp
  ad(1,1)=0.01_dp; ad(2,2)=0.01_dp; v2=0.0_dp; v2(1,1)=3.0_dp; v2(2,2)=3.0_dp
  ivr=riv_gibbs(y,x,z,w,mbg,abg,md,ad,v2,5.0_dp,20,2)
  if (.not.all(ieee_is_finite(ivr%betadraw)) .or. .not.all(ieee_is_finite(ivr%sigmadraw))) &
    error stop "test_iv_dp: IV"
  yd(:,1)=[-2.2_dp,-2.0_dp,-1.8_dp,-2.1_dp,-1.9_dp,-2.05_dp,1.8_dp,2.0_dp,2.2_dp,1.9_dp,2.1_dp,2.05_dp]
  m0=0.0_dp; s0=2.0_dp
  dpr=rdp_gibbs(yd,m0,1.0_dp,5.0_dp,s0,1.0_dp,20,2)
  if (.not.all(ieee_is_finite(dpr%alphadraw)) .or. any(dpr%ncompdraw<1)) error stop "test_iv_dp: DP"
  print *, "test_iv_dp: PASS"
end program test_iv_dp
