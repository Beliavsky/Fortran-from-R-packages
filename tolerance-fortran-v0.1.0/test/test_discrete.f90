program test_discrete
  use tolerance
  implicit none
  type(discrete_tolerance_interval)::di
  type(acceptance_plan)::ap
  integer::fail
  fail=0
  di=bintol_int(12,20,30,0.05_dp,0.9_dp,1,'CP')
  if(di%lower<0 .or. di%upper>30 .or. di%upper<di%lower)call bad('binomial')
  di=poistol_int(25,10,15,0.05_dp,0.9_dp,1,'SC')
  if(di%upper<di%lower)call bad('poisson')
  di=negbintol_int(20,10,15,0.05_dp,0.9_dp,1,'LS')
  if(di%upper<di%lower)call bad('negative binomial')
  ap=acceptance_sampling(50,1000,0.05_dp,0.95_dp,0.01_dp,0.05_dp)
  if(ap%sample_size/=50 .or. ap%acceptance_limit<0)call bad('acceptance')
  if(fail==0)then;print '(a)','test_discrete: PASS';else;error stop 1;end if
contains
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
