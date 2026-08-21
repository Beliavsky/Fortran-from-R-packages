program test_geostat
  use compositions
  implicit none
  real(dp) :: comp(3,3),loc(3,1),bins(2,2),gobs(3,3,3,3),gnew(3,1,3,3),dist
  real(dp), allocatable :: pred(:,:)
  type(variogram_result) :: vr
  integer :: i,j,k
  comp=reshape([0.6_dp,0.3_dp,0.1_dp, 0.5_dp,0.35_dp,0.15_dp, 0.4_dp,0.4_dp,0.2_dp],[3,3],order=[2,1])
  loc(:,1)=[0.0_dp,1.0_dp,2.0_dp]; bins=reshape([0.0_dp,1.1_dp,1.1_dp,2.1_dp],[2,2],order=[2,1])
  vr=logratio_variogram(comp,loc,bins)
  if(any(vr%count<0).or.sum(vr%count)==0) error stop 'variogram counts'
  if(abs(vgram_spherical(0.0_dp,0.1_dp,1.0_dp,2.0_dp))>1e-15_dp) error stop 'spherical origin'
  if(vgram_exponential(2.0_dp,0.1_dp,1.0_dp,2.0_dp)<=0.0_dp) error stop 'exponential variogram'
  gobs=0.0_dp; gnew=0.0_dp
  do i=1,3; do j=1,3
    dist=abs(loc(i,1)-loc(j,1))
    do k=1,3; gobs(i,j,k,k)=dist; end do
  end do; end do
  do i=1,3
    dist=abs(loc(i,1)-loc(1,1)); do k=1,3; gnew(i,1,k,k)=dist; end do
  end do
  pred=compositional_ordinary_kriging(comp,gobs,gnew)
  if(abs(sum(pred(1,:))-1.0_dp)>1e-10_dp) error stop 'kriging closure'
  if(maxval(abs(pred(1,:)-comp(1,:)))>1e-7_dp) then; print *,pred(1,:),comp(1,:); error stop 'kriging interpolation'; end if
  print *, 'test_geostat: PASS'
end program
