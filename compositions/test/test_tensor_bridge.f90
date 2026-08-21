program test_tensor_bridge
  use compositions
  implicit none
  type(tensor_t) :: a,b,c,m
  real(dp) :: ad(12),bd(6)
  integer :: i
  logical :: has_i,has_j,has_k,has_b
  do i=1,12; ad(i)=real(i,dp); end do
  bd=[1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp]
  a=tensor(ad,[2,3,2],['i','j','b'])
  b=tensor(bd,[3,2],['j','k'])
  c=einstein_pair(a,b)
  if(c%rank()/=3) error stop 'named contraction rank'
  if(any(c%shape/=[2,2,2])) error stop 'named contraction shape'
  has_i=.false.; has_j=.false.; has_k=.false.; has_b=.false.
  do i=1,c%rank()
    if(trim(c%axis(i))=='i') has_i=.true.
    if(trim(c%axis(i))=='j') has_j=.true.
    if(trim(c%axis(i))=='k') has_k=.true.
    if(trim(c%axis(i))=='b') has_b=.true.
  end do
  if(.not.has_i.or..not.has_k.or..not.has_b.or.has_j) error stop 'named contraction axes'
  m=mean_tensor(a,['i'])
  if(m%rank()/=2.or.any(m%shape/=[3,2])) error stop 'high-rank named mean'
  if(a%axis_pos('j')/=2) error stop 'dynamic axis lookup'
  print *, 'test_tensor_bridge: PASS'
end program test_tensor_bridge
