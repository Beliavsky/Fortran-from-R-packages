! SPDX-License-Identifier: MIT
program test_psis
  use bayesianou
  implicit none
  real(dp)::ll(100,4)
  type(loo_result)::r
  integer::i,j
  do j=1,4;do i=1,100;ll(i,j)=-0.5_dp*(0.01_dp*i/100.0_dp+j*0.1_dp);end do;end do
  call psis_loo(ll,r)
  if(.not.(r%elpd_loo<0.0_dp))error stop 'bad ELPD'
  if(any(r%pareto_k>5.0_dp))error stop 'bad k'
  print *, 'test_psis: PASS'
end program test_psis
