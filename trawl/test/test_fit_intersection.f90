program test_fit_intersection
  use trawl, only : dp,trawl_fit_result,poisson_fit_result,nb_fit_result,trawl_spec, &
    fit_exptrawl,fit_supigtrawl,fit_lmtrawl,fit_dexptrawl,fit_marginal_poisson, &
    fit_marginal_nb,fit_trawl_intersection,fit_trawl_intersection_exp, &
    set_trawl_seed
  implicit none
  real(dp)::x(8),r,lm1,lm2
  type(trawl_fit_result)::fr
  type(poisson_fit_result)::pr
  type(nb_fit_result)::nr
  type(trawl_spec)::s1,s2
  integer::fail
  fail=0
  x=[1.0_dp,2.0_dp,4.0_dp,3.0_dp,5.0_dp,4.0_dp,6.0_dp,5.0_dp]
  pr=fit_marginal_poisson(x,2.0_dp)
  call chk(pr%v,sum(x)/8.0_dp/2.0_dp,1.0e-13_dp,'poisson fit')
  nr=fit_marginal_nb(x,2.0_dp)
  call assert_true(nr%status==0,'NB marginal fit')
  fr=fit_exptrawl(x)
  call assert_true(fr%status==0 .and. fr%lambda1>0.0_dp,'Exp fit')

  lm1=1.0_dp/2.0_dp;lm2=1.0_dp/0.5_dp
  r=fit_trawl_intersection_exp(2.0_dp,0.5_dp,lm1,lm2)
  call chk(r,min(lm1,lm2),1.0e-10_dp,'Exp intersection')
  s1%kind='DExp';s1%w=0.25_dp;s1%lambda1=0.2_dp;s1%lambda2=2.0_dp
  s2%kind='Exp';s2%lambda1=0.5_dp
  lm1=s1%w/s1%lambda1+(1.0_dp-s1%w)/s1%lambda2;lm2=1.0_dp/s2%lambda1
  r=fit_trawl_intersection(s1,s2,lm1,lm2)
  call assert_true(r>=0.0_dp .and. r<=min(lm1,lm2)+1.0e-10_dp,'generic intersection')

  call set_trawl_seed(2468)
  fr=fit_supigtrawl(x,gmm_lag=3,itermax=50)
  call assert_true(fr%status==0 .and. fr%gamma>0.0_dp,'supIG optimizer')
  fr=fit_lmtrawl(x,gmm_lag=3,itermax=50)
  call assert_true(fr%status==0 .and. fr%alpha>0.0_dp,'LM optimizer')
  fr=fit_dexptrawl(x,gmm_lag=3,itermax=50)
  call assert_true(fr%status==0 .and. fr%lambda1>0.0_dp .and. fr%lambda2>0.0_dp,'DExp optimizer')
  if(fail==0)then;print '(a)','test_fit_intersection: PASS';else;error stop 1;end if
contains
  subroutine chk(a,b,tol,name);real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::name
    if(abs(a-b)>tol)then;print *, 'FAIL ',trim(name),a,b;fail=fail+1;end if;end subroutine
  subroutine assert_true(ok,name);logical,intent(in)::ok;character(len=*),intent(in)::name
    if(.not.ok)then;print *, 'FAIL ',trim(name);fail=fail+1;end if;end subroutine
end program
