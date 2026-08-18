! SPDX-License-Identifier: GPL-3.0-or-later
program test_simulation
  use degreenet_kinds, only : dp
  use degreenet_rng, only : seed_rng
  use degreenet_models, only : MODEL_DP
  use degreenet_simulation, only : sample_model
  implicit none
  integer::x(2000)
  real(dp)::m
  call seed_rng(12345)
  call sample_model(MODEL_DP,[3.5_dp],size(x),x,1,1000)
  if(any(x<1).or.any(x>1000))error stop 1
  m=sum(real(x,dp))/real(size(x),dp)
  if(m<1.0_dp.or.m>2.0_dp)then;print *,'FAIL simulated mean',m;error stop 1;end if
  print *, 'test_simulation: PASS'
end program
