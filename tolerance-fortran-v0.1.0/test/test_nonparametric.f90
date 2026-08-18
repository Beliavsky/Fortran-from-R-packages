program test_nonparametric
  use tolerance
  implicit none
  real(dp)::x(10)=[1._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp,8._dp,9._dp,10._dp]
  type(tolerance_interval)::ti,os,fos
  real(dp)::lo2(2),up2(2)
  integer::fail,n,nopt
  logical::has_fos
  fail=0
  ti=nptol_int(x,0.1_dp,0.8_dp,1,'WILKS');if(ti%upper<ti%lower)call bad('wilks')
  ti=nptol_int(x,0.1_dp,0.8_dp,2,'WALD');if(ti%upper<ti%lower)call bad('wald')
  ti=npbeta_tol_int(x,0.8_dp,2);if(ti%upper<ti%lower)call bad('npbeta')
  call nptol_hm_options(x,0.1_dp,0.8_dp,lo2,up2,nopt)
  if(nopt<1 .or. nopt>2 .or. any(up2(1:nopt)<lo2(1:nopt)))call bad('hm options')
  call nptol_ym_options(x,0.1_dp,0.8_dp,1,os,fos,has_fos)
  if(.not.has_fos .or. os%upper<os%lower .or. fos%upper<fos%lower)call bad('ym options')
  n=distfree_sample_size(0.05_dp,0.95_dp,1);if(n/=59)call bad('distfree n')
  if(np_order(1,0.05_dp,0.95_dp)<1)call bad('np order')
  if(fail==0)then;print '(a)','test_nonparametric: PASS';else;error stop 1;end if
contains
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
