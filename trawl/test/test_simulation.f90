program test_simulation
  use trawl, only : dp,trawl_spec,set_trawl_seed,sim_univariate_trawl,sim_bivariate_trawl
  implicit none
  integer,allocatable::x(:),y(:,:)
  type(trawl_spec)::s1,s2
  integer::st,fail
  real(dp)::m
  fail=0;s1%kind='Exp';s1%lambda1=1.0_dp;s2%kind='Exp';s2%lambda1=0.7_dp
  call set_trawl_seed(777)
  call sim_univariate_trawl(2000.0_dp,x,delta_t=1.0_dp,burnin=20.0_dp, &
    marginal='Poi',spec=s1,v=3.0_dp,status=st)
  call assert_true(st==0 .and. size(x)==2000,'univariate shape')
  call assert_true(minval(x)>=0,'univariate nonnegative')
  m=sum(real(x,dp))/real(size(x),dp)
  call assert_true(abs(m-3.0_dp)<0.35_dp,'univariate mean')
  call sim_univariate_trawl(500.0_dp,x,delta_t=1.0_dp,burnin=10.0_dp, &
    marginal='NegBin',spec=s1,m=2.0_dp,theta=0.4_dp,status=st)
  call assert_true(st==0 .and. minval(x)>=0,'NB trawl simulation')

  call sim_bivariate_trawl(500.0_dp,y,delta_t=1.0_dp,burnin=10.0_dp,marginal='Poi', &
    dependencetype='fullydep',spec1=s1,spec2=s2,v12=2.0_dp,status=st)
  call assert_true(st==0 .and. size(y,1)==500 .and. size(y,2)==2,'bivariate fullydep shape')
  call assert_true(minval(y)>=0,'bivariate fullydep nonnegative')
  call sim_bivariate_trawl(300.0_dp,y,delta_t=1.0_dp,burnin=10.0_dp,marginal='Poi', &
    dependencetype='dep',spec1=s1,spec2=s2,v1=1.0_dp,v2=1.5_dp,v12=0.8_dp,status=st)
  call assert_true(st==0 .and. minval(y)>=0,'bivariate dep simulation')
  if(fail==0)then;print '(a)','test_simulation: PASS';else;error stop 1;end if
contains
  subroutine assert_true(ok,name);logical,intent(in)::ok;character(len=*),intent(in)::name
    if(.not.ok)then;print *, 'FAIL ',trim(name);fail=fail+1;end if;end subroutine
end program
