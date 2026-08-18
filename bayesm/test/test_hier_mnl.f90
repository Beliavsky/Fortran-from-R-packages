program test_hier_mnl
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  type(mnl_data), allocatable :: d(:)
  type(hier_mnl_result) :: fm,dpout
  real(dp) :: z(3,1),dbar(1),ad(1,1),mubar(1),v(1,1),adir(2)
  integer :: i

  call rng_seed(4412)
  allocate(d(3))
  z=1.0_dp
  dbar=0.0_dp
  ad=0.2_dp
  mubar=0.0_dp
  v=4.0_dp
  adir=1.0_dp

  do i=1,3
    allocate(d(i)%y(4),d(i)%x(8,1))
    d(i)%y=[1,2,1,2]
    d(i)%x(:,1)=[0.0_dp,1.0_dp, 0.2_dp,1.2_dp, -0.2_dp,0.8_dp, 0.1_dp,1.1_dp]
  end do

  fm=rhier_mnl_rw_mixture(d,z,dbar,ad,mubar,0.2_dp,5.0_dp,v,adir,0.15_dp,10,2,.true.)
  if (.not.all(ieee_is_finite(fm%betadraw))) error stop "test_hier_mnl: finite mix beta"
  if (.not.all(ieee_is_finite(fm%llike))) error stop "test_hier_mnl: finite mix ll"
  if (maxval(abs(sum(fm%probdraw,dim=2)-1.0_dp))>1.0e-10_dp) error stop "test_hier_mnl: finite mix p"

  dpout=rhier_mnl_dp(d,z,dbar,ad,mubar,0.2_dp,5.0_dp,v,1.0_dp,3,0.15_dp,10,2)
  if (.not.all(ieee_is_finite(dpout%betadraw))) error stop "test_hier_mnl: dp beta"
  if (maxval(abs(sum(dpout%probdraw,dim=2)-1.0_dp))>1.0e-10_dp) error stop "test_hier_mnl: dp p"

  print *, "test_hier_mnl: PASS"
end program test_hier_mnl
