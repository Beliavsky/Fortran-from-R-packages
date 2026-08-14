program test_conservative
  use anmc
  implicit none
  integer,parameter::n=5
  real(dp)::pn(n),mu(n),sigma(n,n),e(n,1),expected_joint
  integer::i
  type(conservative_result)::r
  pn=[0.99_dp,0.98_dp,0.97_dp,0.80_dp,0.60_dp]
  sigma=0.0_dp
  do i=1,n
    sigma(i,i)=1.0_dp
    mu(i)=normal_quantile(pn(i))
    e(i,1)=real(i,dp)
  end do
  r=conservative_estimate(0.90_dp,mu,sigma,e,0.0_dp,pn=pn,excursion_type='>', &
                          prob_control=genz_bretz(maxpts=100000,batches=16,abseps=1e-6_dp))
  if(.not.r%ok) error stop 'conservative estimate failed'
  if(count(r%set)/=3) then
    write(*,*) 'unexpected set size ',count(r%set)
    error stop 1
  end if
  if(.not.all(r%set(1:3)) .or. any(r%set(4:5))) error stop 'wrong conservative set membership'
  expected_joint=product(pn(1:3))
  if(abs(r%probability-expected_joint)>1.0e-3_dp) then
    write(*,*) 'probability mismatch ',r%probability,expected_joint
    error stop 1
  end if
  if(abs(r%level-0.80_dp)>1.0e-12_dp) error stop 'wrong conservative level'
  print *, 'test_conservative: PASS'
end program test_conservative
