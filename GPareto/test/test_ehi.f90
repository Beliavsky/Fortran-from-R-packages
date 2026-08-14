program test_ehi
  use gpareto, only : dp, ehi_2d_values, hypervolume_improvement
  implicit none
  real(dp) :: front(3,2),ref(2),mu(2,2),sd(2,2),det
  real(dp),allocatable :: e(:)
  front=reshape([1.0_dp,2.0_dp,4.0_dp, 4.0_dp,2.0_dp,1.0_dp],[3,2]);ref=[5.0_dp,5.0_dp]
  mu(1,:)=[2.5_dp,2.5_dp];sd(1,:)=[0.5_dp,0.7_dp]
  mu(2,:)=[1.5_dp,1.5_dp];sd(2,:)=0.0_dp
  call ehi_2d_values(front,ref,mu,sd,e)
  if(abs(e(1)-0.209242250192_dp)>2.0e-6_dp) then
    print *, 'got EHI=',e(1)
    error stop 'analytical EHI regression failed'
  end if
  det=hypervolume_improvement(mu(2,:),front,ref)
  if(abs(e(2)-det)>1.0e-12_dp) error stop 'deterministic EHI failed'
  print *, 'test_ehi PASS'
end program test_ehi
