program test_distribution
  use dirichletreg, only : dp, ddirichlet_log, rdirichlet, seed_rng, prepare_composition, to_ternary
  implicit none
  real(dp) :: x(3,3), alpha(3), amat(3,3), ld(3), ref(3), sim(20000,3), means(3)
  real(dp) :: raw(2,3), prep(2,3), xy(2,2)
  logical :: normed, transformed
  integer :: stat, failures

  failures=0
  x=reshape([0.2_dp,0.1_dp,0.4_dp, 0.3_dp,0.7_dp,0.4_dp, 0.5_dp,0.2_dp,0.2_dp],[3,3])
  alpha=[2.0_dp,3.0_dp,4.0_dp]
  ref=[2.022871190191443_dp,0.275447534783440_dp,0.542510320032486_dp]
  call ddirichlet_log(x,alpha,ld,stat)
  if(stat/=0 .or. maxval(abs(ld-ref))>2.0e-12_dp) failures=failures+1

  amat=reshape([2.0_dp,1.2_dp,4.0_dp, 3.0_dp,2.3_dp,1.5_dp, 4.0_dp,3.4_dp,2.0_dp],[3,3])
  ref=[2.022871190191443_dp,0.444755596085537_dp,1.046931531172281_dp]
  call ddirichlet_log(x,amat,ld,stat)
  if(stat/=0 .or. maxval(abs(ld-ref))>2.0e-12_dp) failures=failures+1

  call seed_rng(12345)
  call rdirichlet(size(sim,1),alpha,sim,stat)
  means=sum(sim,dim=1)/real(size(sim,1),dp)
  if(stat/=0 .or. maxval(abs(sum(sim,dim=2)-1.0_dp))>2.0e-14_dp) failures=failures+1
  if(maxval(abs(means-alpha/sum(alpha)))>0.01_dp) failures=failures+1

  raw=reshape([2.0_dp,0.0_dp, 3.0_dp,1.0_dp, 5.0_dp,1.0_dp],[2,3])
  call prepare_composition(raw,prep,normalized=normed,transformed=transformed,stat=stat)
  if(stat/=0 .or. .not.normed .or. .not.transformed) failures=failures+1
  if(maxval(abs(sum(prep,dim=2)-1.0_dp))>1.0e-14_dp .or. any(prep<=0.0_dp)) failures=failures+1

  call to_ternary(prep,xy,stat)
  if(stat/=0) failures=failures+1

  if(failures==0) then
    print '(a)', 'test_distribution: PASS'
  else
    print '(a,i0)', 'test_distribution: FAIL ', failures
    error stop 1
  end if
end program test_distribution
