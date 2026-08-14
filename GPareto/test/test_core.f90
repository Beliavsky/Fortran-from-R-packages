program test_core
  use gpareto, only : dp, hypervolume, hypervolume_improvement, nondominated_points, nondom_set
  use gpareto, only : zdt1, p1_test
  implicit none
  real(dp) :: p(4,2), ref(2), hv, imp, p3(4,3), ref3(3)
  real(dp), allocatable :: front(:,:), y(:,:)
  logical, allocatable :: keep(:)
  p = reshape([1.0_dp,2.0_dp,4.0_dp,4.5_dp, 4.0_dp,2.0_dp,1.0_dp,4.5_dp],[4,2])
  ref=[5.0_dp,5.0_dp]
  call nondominated_points(p,front)
  if(size(front,1)/=3) error stop 'nondominated_points failed'
  hv=hypervolume(front,ref)
  if(abs(hv-11.0_dp)>1.0e-12_dp) error stop 'hypervolume 2D failed'
  imp=hypervolume_improvement([1.0_dp,1.0_dp],front,ref)
  if(abs(imp-5.0_dp)>1.0e-12_dp) error stop 'hypervolume improvement failed'
  p3=reshape([1.0_dp,2.0_dp,3.0_dp,1.5_dp, 3.0_dp,2.0_dp,1.0_dp,2.5_dp, &
    4.0_dp,3.0_dp,2.0_dp,2.5_dp],[4,3])
  ref3=[5.0_dp,5.0_dp,5.0_dp]
  if(abs(hypervolume(p3,ref3)-35.375_dp)>1.0e-12_dp) error stop 'hypervolume 3D failed'
  call nondom_set(reshape([3.0_dp,0.5_dp,3.0_dp,0.5_dp],[2,2]),front,keep)
  if(keep(1)) error stop 'nondom_set dominated point failed'
  call zdt1(reshape([0.25_dp,0.0_dp],[1,2]),y)
  if(abs(y(1,1)-0.25_dp)>1.0e-12_dp .or. abs(y(1,2)-0.5_dp)>1.0e-12_dp) error stop 'ZDT1 failed'
  call p1_test(reshape([0.5_dp,0.5_dp],[1,2]),y)
  if(any(.not.(abs(y)<huge(1.0_dp)))) error stop 'P1 failed'
  print *, 'test_core PASS'
end program test_core
