program test_hierarchical
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  type(reg_data), allocatable :: d(:)
  type(hier_linear_result) :: hr,hm
  real(dp) :: z(3,1),db(1,2),a(1,1),v(2,2),ssq(3),dvec(2),ad(2,2),mubar(2),adir(2)
  integer :: i,j
  call rng_seed(991)
  allocate(d(3)); z=1.0_dp; db=0.0_dp; a=0.1_dp; v=0.0_dp; v(1,1)=4.0_dp; v(2,2)=4.0_dp; ssq=1.0_dp
  do i=1,3
    allocate(d(i)%y(6),d(i)%x(6,2)); d(i)%x(:,1)=1.0_dp
    do j=1,6; d(i)%x(j,2)=real(j-3,dp); end do
    d(i)%y=0.3_dp*real(i,dp)+0.4_dp*d(i)%x(:,2)
  end do
  hr=rhier_linear_model(d,z,db,a,6.0_dp,v,5.0_dp,ssq,20,2)
  if (.not.all(ieee_is_finite(hr%betadraw))) error stop "test_hierarchical: linear"
  dvec=0.0_dp; ad=0.0_dp; ad(1,1)=0.1_dp; ad(2,2)=0.1_dp; mubar=0.0_dp; adir=1.0_dp
  hm=rhier_linear_mixture(d,z,dvec,ad,mubar,0.1_dp,6.0_dp,v,5.0_dp,ssq,adir,10,2,.true.)
  if (.not.all(ieee_is_finite(hm%betadraw))) error stop "test_hierarchical: mixture"
  if (maxval(abs(sum(hm%probdraw,dim=2)-1.0_dp))>1.0e-12_dp) error stop "test_hierarchical: mix p"
  print *, "test_hierarchical: PASS"
end program test_hierarchical
