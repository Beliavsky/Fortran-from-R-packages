program test_random
  use tmvtnorm
  implicit none
  real(dp)::mu(2),s(2,2),lo(2),up(2),d(3,2),a(3),b(3),start(2),h(2,2)
  real(dp),allocatable::x(:,:),y(:,:),z(:,:),ts(:,:)
  integer::rowind(4),colptr(3)
  real(dp)::val(4)
  mu=0.0_dp
  s=reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
  lo=-1.0_dp
  up=1.0_dp
  x=rtmvnorm_gibbs(500,mu,s,lo,up,burnin=50,seed=123)
  if(any(x<spread(lo,1,500)) .or. any(x>spread(up,1,500))) error stop 'gibbs support'
  if(maxval(abs(sum(x,dim=1)/500.0_dp))>0.13_dp) error stop 'gibbs mean'
  y=rtmvnorm_rejection(100,mu,s,lo,up,seed=321)
  if(any(y<spread(lo,1,100)) .or. any(y>spread(up,1,100))) error stop 'rejection support'

  d=reshape([1.0_dp,1.0_dp,0.5_dp,1.0_dp,-1.0_dp,-1.0_dp],[3,2])
  a=0.0_dp
  b=1.0_dp
  start=[0.4_dp,0.2_dp]
  z=rtmvnorm2(150,mu,s,d,a,b,algorithm=algorithm_gibbs,burnin=20,start=start,seed=99)
  if(any(matmul(z,transpose(d))<spread(a,1,150)-1e-10_dp) .or. &
     any(matmul(z,transpose(d))>spread(b,1,150)+1e-10_dp)) error stop 'linear constraints'

  h=reshape([1.19047619047619_dp,-0.476190476190476_dp,-0.476190476190476_dp,1.19047619047619_dp],[2,2])
  rowind=[1,2,1,2]
  colptr=[1,3,5]
  val=[h(1,1),h(2,1),h(1,2),h(2,2)]
  z=rtmvnorm_sparse_csc(100,mu,rowind,colptr,val,lo,up,burnin=10,seed=77)
  if(any(z<spread(lo,1,100)) .or. any(z>spread(up,1,100))) error stop 'sparse support'

  ts=rtmvt_gibbs(100,mu,s,5.0_dp,lo,up,burnin=10,seed=41)
  if(any(ts<spread(lo,1,100)) .or. any(ts>spread(up,1,100))) error stop 't gibbs support'
  print *, 'test_random: ok'
end program
