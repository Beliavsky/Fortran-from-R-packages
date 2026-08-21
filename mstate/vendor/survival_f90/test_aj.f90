program test_aj
  use survival
  implicit none
  real(dp)::time(4)
  integer::from(4),to(4)
  real(dp),allocatable::p(:,:),ut(:)
  time=[1._dp,2._dp,3._dp,4._dp]
  from=[1,1,1,2];to=[2,2,3,3]
  call aalen_johansen(time,from,to,3,p,ut)
  if(any(p< -1e-12_dp)) error stop 'aj negative'
  if(abs(sum(p(:,size(p,2)))-1._dp)>1e-10_dp) error stop 'aj sum'
  print *, 'test_aj PASS'
end program
