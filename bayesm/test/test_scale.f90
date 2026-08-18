program test_scale
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  integer :: x(5,2)
  real(dp) :: v(2,2),mubar(2),am(2,2),gs(5),g11(5),g22(5),g12(5),vl(2,2),ge(3)
  type(scale_usage_result) :: out
  call rng_seed(404)
  x=reshape([1,2,2,3,3, 2,2,3,3,4],[5,2])
  v=0.0_dp; v(1,1)=4.0_dp; v(2,2)=4.0_dp; mubar=2.0_dp; am=0.0_dp; am(1,1)=0.05_dp; am(2,2)=0.05_dp
  gs=[0.5_dp,0.8_dp,1.0_dp,1.3_dp,1.8_dp]; g11=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  g22=[0.3_dp,0.5_dp,0.8_dp,1.0_dp,1.5_dp]; g12=[-0.4_dp,-0.2_dp,0.0_dp,0.2_dp,0.4_dp]
  vl=0.0_dp; vl(1,1)=10.0_dp; vl(2,2)=2.0_dp; ge=[-0.05_dp,0.0_dp,0.05_dp]
  out=rscale_usage(x,4,5.0_dp,v,mubar,am,gs,g11,g22,g12,10,vl,ge,4,1,10)
  if (.not.all(ieee_is_finite(out%mudraw)) .or. .not.all(ieee_is_finite(out%sigmadraw))) &
    error stop "test_scale: finite"
  if (any(out%sdraw<=0.0_dp)) error stop "test_scale: sigma"
  print *, "test_scale: PASS"
end program test_scale
