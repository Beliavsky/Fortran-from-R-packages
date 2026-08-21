program test_gof
  use compositions
  implicit none
  real(dp) :: x(3,2),y(3,2),stat,ps(5),rep(20)
  integer :: dat(5)
  x=reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[3,2],order=[2,1])
  y=x
  stat=kernel_similarity_statistic(x,y,0.5_dp)
  if(abs(stat-1.0878648274324869_dp)>1e-12_dp) then
    print *,stat; error stop 'kernel source statistic'
  end if
  ps=exp(-1.0_dp)*[1.0_dp,1.0_dp,0.5_dp,1.0_dp/6.0_dp,1.0_dp/24.0_dp]
  dat=[0,1,1,2,0]
  stat=poisson_ks_statistic(dat,ps)
  if(stat<0.0_dp.or.stat>1.0_dp) error stop 'Poisson KS range'
  call rng_seed(4567); rep=poisson_ks_sample(30,ps,20)
  if(any(rep<0.0_dp).or.any(rep>1.0_dp)) error stop 'Poisson KS sample range'
  print *, 'test_gof: PASS'
end program
