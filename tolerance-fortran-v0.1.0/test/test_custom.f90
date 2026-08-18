program test_custom
  use tolerance
  implicit none
  integer::xd(8)=[0,0,1,1,2,3,4,7], xp(10)=[0,1,0,2,3,1,4,0,2,1]
  integer::xz(12)=[1,1,1,2,2,3,1,4,2,5,1,3]
  real(dp)::th,s,b
  type(discrete_tolerance_interval)::di
  integer::fail
  fail=0
  th=dpareto_mle(xd);if(th<=0._dp .or. th>=1._dp)call bad('dpareto mle')
  di=dparetotol_int(xd,alpha=0.1_dp,p=0.9_dp,theta_hat=th);if(di%lower<0)call bad('dpareto tol')
  th=poislind_mle(xp);if(th<=0._dp)call bad('poislind mle')
  di=zipftol_int(xz,nmax=5,alpha=0.1_dp,p=0.9_dp,dist='Zipf',s_hat=s,b_hat=b)
  if(s<=0._dp .or. di%lower<1 .or. di%upper>5)call bad('zipf')
  if(fail==0)then;print '(a)','test_custom: PASS';else;error stop 1;end if
contains
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
