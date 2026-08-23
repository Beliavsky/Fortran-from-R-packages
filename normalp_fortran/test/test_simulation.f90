program test_simulation
  use normalp
  implicit none
  type(simulation_summary) :: a,b
  real(dp) :: beta(2)
  integer, allocatable :: seed(:)
  integer :: nseed
  call random_seed(size=nseed); allocate(seed(nseed)); seed=24680; call random_seed(put=seed)
  call simul_mp(80,20,a,mu=1.0_dp,sigmap=0.8_dp,p=2.0_dp)
  if(size(a%draws,1)/=20 .or. size(a%draws,2)/=5) error stop 1
  if(abs(a%means(1)-1.0_dp)>0.2_dp) error stop 2
  beta=[0.5_dp,-1.2_dp]
  call simul_lmp(80,10,beta,b,sigmap=0.2_dp,p=2.0_dp,estimate_shape=.false.)
  if(abs(b%means(1)-beta(1))>0.15_dp .or. abs(b%means(2)-beta(2))>0.2_dp) error stop 3
  print '(a)', 'test_simulation: PASS'
end program
