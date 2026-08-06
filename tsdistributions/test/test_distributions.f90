program test_distributions
  use tsdistributions
  implicit none
  character(len=8), parameter :: names(10)=[character(len=8) :: 'norm','std','snorm','sstd','ged','sged','nig','gh','jsu','ghst']
  type(distribution_parameters)::p
  real(dp)::q,prob,d
  integer::i
  p=distribution_parameters(mu=0.1_dp,sigma=1.2_dp,skew=0.2_dp,shape=8.2_dp,lambda=-0.5_dp)
  do i=1,size(names)
    select case(trim(names(i)))
    case('std');p%skew=1.0_dp;p%shape=8.0_dp
    case('snorm');p%skew=1.3_dp
    case('sstd');p%skew=1.3_dp;p%shape=8.0_dp
    case('ged');p%skew=1.0_dp;p%shape=1.5_dp
    case('sged');p%skew=1.2_dp;p%shape=1.5_dp
    case('nig');p%skew=0.2_dp;p%shape=1.5_dp
    case('gh');p%skew=0.2_dp;p%shape=2.0_dp;p%lambda=-0.5_dp
    case('jsu');p%skew=0.2_dp;p%shape=1.2_dp
    case('ghst');p%skew=0.2_dp;p%shape=8.2_dp
    end select
    prob=0.3_dp;q=qdist(trim(names(i)),prob,p);d=ddist(trim(names(i)),q,p)
    call assert_true(abs(pdist(trim(names(i)),q,p)-prob)<2.0e-6_dp,'CDF/quantile inversion '//trim(names(i)))
    call assert_true(d>0.0_dp,'positive density '//trim(names(i)))
    call assert_true(abs(exp(ddist(trim(names(i)),q,p,.true.))-d)<1.0e-10_dp,'log density '//trim(names(i)))
  end do
  call assert_true(abs(pdist('norm',0.0_dp,distribution_parameters())-0.5_dp)<1.0e-12_dp,'normal CDF')
  q=qghyp(0.4_dp,alpha=1.4_dp,beta_value=0.2_dp,delta=1.1_dp,mu=0.1_dp,lambda=-0.5_dp)
  call assert_true(abs(pghyp(q,alpha=1.4_dp,beta_value=0.2_dp,delta=1.1_dp,mu=0.1_dp,lambda=-0.5_dp)-0.4_dp)<2.0e-6_dp,'raw GH inversion')
  print '(a)','test_distributions: PASS'
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine
end program
