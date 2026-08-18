program test_continuous
  use tolerance
  implicit none
  real(dp) :: x(5)=[1._dp,2._dp,3._dp,4._dp,5._dp], sh,sc
  type(tolerance_interval)::ti
  integer::fail
  fail=0
  ti=exptol_int(x,0.05_dp,0.9_dp,1,.false.)
  call chk(ti%lower,0.17265575460850247_dp,2e-10_dp,'exp lower')
  call chk(ti%upper,17.53104279738982_dp,2e-9_dp,'exp upper')
  ti=uniftol_int(x,0.05_dp,0.9_dp,2)
  if(.not.(ti%lower<ti%upper))call bad('uniform ordering')
  ti=laptol_int(x,0.05_dp,0.9_dp,1)
  if(.not.(ti%lower<3._dp .and. ti%upper>3._dp))call bad('laplace center')
  ti=gamtol_int(x,0.05_dp,0.9_dp,1,'HE',50,.false.,sh,sc)
  if(.not.(sh>0._dp .and. sc>0._dp .and. ti%upper>ti%lower))call bad('gamma fit')
  ti=cautol_int(x,0.05_dp,0.9_dp,1)
  if(ti%upper<=ti%lower)call bad('cauchy')
  if(fail==0)then;print '(a)','test_continuous: PASS';else;error stop 1;end if
contains
  subroutine chk(a,b,t,nm)
    real(dp),intent(in)::a,b,t;character(len=*),intent(in)::nm
    if(abs(a-b)>t)then;print *,trim(nm),a,b;fail=fail+1;end if
  end subroutine
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
