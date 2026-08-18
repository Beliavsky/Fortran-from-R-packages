program test_multivariate
  use tolerance
  implicit none
  real(dp)::x(8,2),k,depth(8),lims(2,2)
  integer::tv(2)
  type(mv_tolerance_region)::r
  integer::i,fail
  fail=0
  do i=1,8;x(i,:)=[real(i,dp),0.5_dp*real(i,dp)+real(mod(i,2),dp)];depth(i)=1._dp/real(1+abs(i-4),dp);end do
  k=mvtol_factor(20,2,0.05_dp,0.9_dp,'AM');if(abs(k-7.382647443824329_dp)>2e-7_dp)call bad('AM factor')
  r=mvtol_region(x,0.1_dp,0.8_dp,'AM');if(r%radius2<=0._dp)call bad('region')
  call npmvtol_region(x,depth,lims,alpha=0.1_dp,content=0.7_dp);if(any(lims(:,2)<lims(:,1)))call bad('npmv')
  tv=[1,0];call npmvtol_region(x,depth,lims,alpha=0.1_dp,content=0.7_dp,typevec=tv,lower_inf=-999._dp)
  if(lims(1,1)/=-999._dp)call bad('npmv semispace')
  if(fail==0)then;print '(a)','test_multivariate: PASS';else;error stop 1;end if
contains
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
