program test_zero_counts
  use compositions
  implicit none
  real(dp) :: x(3),d(3),y(3),p(3,3),sp(3,3),pr(3)
  logical :: has(3),hm(2,3)
  integer :: cnt(1000,3),tot(1000)
  x=[0.0_dp,0.3_dp,0.7_dp]; d=[0.03_dp,0.03_dp,0.03_dp]
  y=zero_replace(x,d)
  if(abs(y(1)-0.02_dp)>1e-12_dp) error stop 'zero replace'
  has=[.true.,.false.,.true.]; p=missing_projector_acomp(has)
  if(maxval(abs(matmul(p,p)-p))>1e-12_dp) error stop 'projector idempotence'
  hm=reshape([.true.,.true.,.true.,.false.,.false.,.true.],[2,3],order=[2,1])
  sp=sum_missing_projector_acomp(hm)
  if(maxval(abs(sp-transpose(sp)))>1e-12_dp) error stop 'sum projector symmetry'
  call rng_seed(9876)
  cnt=rmultinom_composition(1000,[0.2_dp,0.3_dp,0.5_dp],100)
  tot=count_totals(cnt); if(any(tot/=100)) error stop 'multinomial total'
  pr=sum(count_proportions(cnt),dim=1)/1000.0_dp
  if(maxval(abs(pr-[0.2_dp,0.3_dp,0.5_dp]))>0.02_dp) error stop 'multinomial mean'
  print *, 'test_zero_counts: PASS'
end program
