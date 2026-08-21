program test_missing_geostat
  use compositions
  implicit none
  real(dp) :: comp(3,3),zmiss(3,3),loc(3),gobs(3,3,3,3),gnew(3,1,3,3),f(3,1),fn(1,1),dist
  real(dp), allocatable :: p0(:,:),p1(:,:),kv(:,:,:)
  integer :: i,j,k
  comp=reshape([0.6_dp,0.3_dp,0.1_dp, 0.5_dp,0.35_dp,0.15_dp, 0.4_dp,0.4_dp,0.2_dp],[3,3],order=[2,1])
  loc=[0.0_dp,1.0_dp,2.0_dp]; gobs=0.0_dp; gnew=0.0_dp; f=1.0_dp; fn=1.0_dp
  do i=1,3; do j=1,3
    dist=abs(loc(i)-loc(j)); do k=1,3; gobs(i,j,k,k)=dist; end do
  end do; end do
  do i=1,3
    dist=abs(loc(i)-loc(1)); do k=1,3; gnew(i,1,k,k)=dist; end do
  end do
  p0=compositional_ordinary_kriging(comp,gobs,gnew)
  p1=compositional_general_kriging(comp,f,gobs,fn,gnew,kv)
  if(maxval(abs(p0-p1))>1.0e-9_dp) error stop 'general kriging complete parity'
  if(maxval(abs(kv(1,:,:)))>1.0e-8_dp) error stop 'general kriging interpolation covariance'
  zmiss=comp; zmiss(2,1)=0.0_dp
  p1=compositional_general_kriging(zmiss,f,gobs,fn,gnew)
  if(abs(sum(p1(1,:))-1.0_dp)>1.0e-12_dp.or.any(p1<=0.0_dp)) error stop 'missing kriging closure'
  if(maxval(abs(p1(1,:)-comp(1,:)))>1.0e-7_dp) error stop 'missing kriging interpolation'
  print *, 'test_missing_geostat: PASS'
end program
