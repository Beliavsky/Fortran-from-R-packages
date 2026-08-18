program test_fiducial
  use tolerance
  implicit none
  type(tolerance_interval)::ti
  real(dp)::zci(2),zpi,zti,zta,lci(2),lpi,lti,lta
  real(dp)::x(10)=[0._dp,0._dp,1._dp,1.2_dp,0.7_dp,2._dp,1.4_dp,0.9_dp,1.8_dp,0.5_dp]
  integer::fail
  fail=0
  ti=fidbintol_int(5,8,10,12,difference,alpha=0.1_dp,p=0.8_dp,kouter=30,binner=40)
  if(ti%upper<ti%lower)call bad('fid bin')
  call semiconttol_int(x,zci,zpi,zti,zta,lci,lpi,lti,lta,0.1_dp,0.8_dp,80)
  if(zci(2)<zci(1) .or. lci(2)<lci(1))call bad('semicont')
  if(fail==0)then;print '(a)','test_fiducial: PASS';else;error stop 1;end if
contains
  real(dp) function difference(a,b) result(v)
    real(dp),intent(in)::a,b;v=a-b
  end function
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
