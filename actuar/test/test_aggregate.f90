! SPDX-License-Identifier: GPL-2.0-or-later
program test_aggregate
  use actuar, only : dp, aggregate_distribution, panjer_poisson, &
    exact_compound, aggregate_var, aggregate_cte
  implicit none
  type(aggregate_distribution) :: panjer,exact
  real(dp) :: severity(3),frequency(7)
  integer :: n

  severity=[0.2_dp,0.5_dp,0.3_dp]
  panjer=panjer_poisson(severity,1.5_dp,tolerance=0.999999_dp,max_terms=80)
  if(.not.panjer%ok) error stop 1
  call assert_close(panjer%probability(1),0.3011942119122021_dp,2.0e-13_dp)
  call assert_close(panjer%probability(2),0.2258956589341516_dp,2.0e-13_dp)
  call assert_close(panjer%probability(3),0.2202482674607978_dp,2.0e-13_dp)

  do n=0,6
    frequency(n+1)=exp(-1.5_dp)*1.5_dp**n/gamma(real(n+1,dp))
  end do
  exact=exact_compound(severity,frequency)
  if(.not.exact%ok) error stop 1
  call assert_close(exact%probability(1),sum(frequency*severity(1)** &
    [(real(n,dp),n=0,6)]),1.0e-12_dp)
  if(aggregate_var(panjer,0.75_dp)<2.0_dp) error stop 1
  if(aggregate_cte(panjer,0.9_dp)<=aggregate_var(panjer,0.9_dp)) error stop 1

  print '(a)', 'test_aggregate: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp), intent(in) :: actual,expected,tol
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected))) then
      print '(a,3es24.15)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_aggregate
