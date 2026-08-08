program test_fitness
  use calibrar, only : dp, fitness_norm2, fitness_pois, fitness_lnorm2, fitness_rangeq
  implicit none
  real(dp) :: obs(3),sim(3),q,p
  obs=[1.0_dp,2.0_dp,3.0_dp];sim=[1.0_dp,1.0_dp,4.0_dp]
  if(abs(fitness_norm2(obs,sim)-2.0_dp)>1.0e-12_dp) error stop "norm2 failed"
  if(.not.(fitness_pois(obs,sim)>0.0_dp)) error stop "pois failed"
  if(.not.(fitness_lnorm2(obs,sim)>0.0_dp)) error stop "lnorm2 failed"
  p=fitness_rangeq(obs,sim,qout=q)
  if(.not.(q>0.0_dp .and. p>=0.0_dp)) error stop "rangeq failed"
  print *, "PASS test_fitness"
end program test_fitness
