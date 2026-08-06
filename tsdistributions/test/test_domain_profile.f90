program test_domain_profile
  use tsdistributions
  implicit none
  type(authorized_domain_result)::domain
  type(parameter_specification)::spec
  type(profile_summary)::profile
  type(rng_state)::rng
  real(dp)::y(10)
  integer::sizes(1),i
  y=[(real(i,dp),i=1,10)]
  domain=authorized_domain('sstd',max_kurt=10.0_dp,n=6)
  call assert_true(domain%status==tsd_success.and.size(domain%skewness)==6,'authorized domain')
  call assert_true(all(domain%kurtosis>=2.9_dp),'domain kurtosis')
  spec=distribution_modelspec(y,'norm')
  spec%parameters%mu=0.0_dp;spec%parameters%sigma=1.0_dp
  call seed_rng(rng,9901_i8);sizes=[80]
  profile=tsprofile(spec,sizes,3,rng,max_iterations=300)
  call assert_true(profile%status==tsd_success.and.profile%successful_fits>0,'simulation profile')
  call assert_true(.not.valid_distribution('banana'),'invalid distribution')
  print '(a)','test_domain_profile: PASS'
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine
end program
