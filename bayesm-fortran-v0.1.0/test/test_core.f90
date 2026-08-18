program test_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  use bayesm_kinds, only: log2pi
  implicit none
  real(dp) :: x2(2),mu2(2),root2(2,2),val
  real(dp), allocatable :: xd(:,:),design(:,:)
  real(dp) :: y(6),xr(6,2),bbar(2),ap(2,2)
  type(unireg_result) :: ur
  integer :: info,i
  call rng_seed(12345)
  x2=0.0_dp; mu2=0.0_dp; root2=0.0_dp; root2(1,1)=1.0_dp; root2(2,2)=1.0_dp
  val=lnd_mvn(x2,mu2,root2)
  if (abs(val+log2pi)>1.0e-13_dp) error stop "test_core: lnd_mvn"
  allocate(xd(2,1)); xd(:,1)=[2.0_dp,3.0_dp]
  call create_x(3,xd=xd,xout=design,info=info)
  if (info/=0 .or. size(design,1)/=6 .or. size(design,2)/=4) error stop "test_core: create_x shape"
  xr(:,1)=1.0_dp
  do i=1,6; xr(i,2)=real(i-1,dp); end do
  y=1.0_dp+2.0_dp*xr(:,2); bbar=0.0_dp; ap=0.0_dp; ap(1,1)=0.01_dp; ap(2,2)=0.01_dp
  ur=runireg(y,xr,bbar,ap,5.0_dp,1.0_dp,20,2)
  if (size(ur%betadraw,1)/=10 .or. .not.all(ieee_is_finite(ur%betadraw))) error stop "test_core: runireg"
  print *, "test_core: PASS"
end program test_core
