program test_em
  use mixsqp
  implicit none
  integer,parameter::n=5,m=3
  real(dp)::L(n,m),w(n),x(m),f0,f1
  L=reshape([1._dp,.2_dp,.7_dp,.4_dp,.9_dp, .2_dp,1._dp,.3_dp,.8_dp,.5_dp, .5_dp,.6_dp,1.1_dp,.2_dp,.7_dp],[n,m])
  w=1._dp/real(n,dp);x=1._dp/real(m,dp)
  f0=mixobjective(L,x,w)
  call mixem_update(L,w,x)
  f1=mixobjective(L,x,w)
  if(f1>f0+1e-12_dp) error stop 'EM increased objective'
  if(abs(sum(x)-1._dp)>1e-12_dp) error stop 'EM weights not normalized'
  print *, 'test_em: PASS'
end program
