program test_data_utils
  use relsurv
  implicit none
  type(split_data_type)::sp
  type(ratetable_type)::tab
  real(dp)::stop(2),cut(1),men(2,2),women(2,2),yc(2),x(1,3),t(1),s(1)
  integer::event(2)
  stop=[3.0_dp,1.5_dp]; cut=2.0_dp; event=[1,0]
  call survsplit_counting(stop=stop,event=event,cut=cut,out=sp)
  if(size(sp%stop)/=3) then; print *, 'FAIL split size'; error stop 1; end if
  if(sum(sp%event)/=1) then; print *, 'FAIL split event'; error stop 1; end if
  men=0.99_dp; women=0.98_dp; yc=[0.0_dp,365.241_dp]
  tab=transrate(men,women,yc)
  x(1,:)=[0.0_dp,1.0_dp,0.0_dp]; t=365.241_dp
  s=expected_survival(tab,x,t)
  if(abs(s(1)-0.99_dp)>2.0e-8_dp) then; print *, 'FAIL transrate',s; error stop 1; end if
  print *, 'test_data_utils: PASS'
end program
