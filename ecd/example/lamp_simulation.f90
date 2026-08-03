! SPDX-License-Identifier: Artistic-2.0
program lamp_simulation
  use ecd_api
  implicit none
  type(rng_state) :: rng
  type(lamp_model) :: model
  type(lamp_result) :: simulation
  type(sample_statistics) :: stats
  integer :: status

  call rng_seed(rng,20260729_i8)
  model=lamp_new(lambda=4.0_dp,beta=0.0_dp,random_walk=22,t_infinity=50.0_dp, &
    random_count=20000,n_lower=0.0_dp,n_upper=1000.0_dp,status=status)
  if(status/=ecd_ok) error stop 'invalid LAMP model'
  call lamp_simulate(model,rng,500,simulation,drop=5)
  if(simulation%status/=ecd_ok .or. size(simulation%z)/=500) error stop 'LAMP simulation failed'
  stats=sample_stats(simulation%z)
  write(*,'(a,i0)') 'observations: ',size(simulation%z)
  write(*,'(a,f14.6)') 'mean:         ',stats%mean
  write(*,'(a,f14.6)') 'standard dev: ',stats%sd
  write(*,'(a,f14.6)') 'minimum:      ',minval(simulation%z)
  write(*,'(a,f14.6)') 'maximum:      ',maxval(simulation%z)
end program lamp_simulation
