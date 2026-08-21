program test_expprep2
  use relsurv
  implicit none
  type(ratetable_type) :: tab
  type(net_summary_type) :: s1,s2
  integer :: dims(1),factor(1),ncuts(1),status(2)
  real(dp) :: cuts(1,1),rates(1),x(2,1),y(2),times(2),t(2),surv(2)
  dims=1; factor=1; ncuts=0; cuts=0.0_dp; rates=0.02_dp
  tab=make_ratetable(dims,factor,cuts,ncuts,rates)
  x=1.0_dp; t=[1.0_dp,2.0_dp]
  surv=expprep2_expected(tab,x,t)
  if(maxval(abs(surv-exp(-0.02_dp*t)))>1.0e-12_dp) error stop 1
  y=[1.0_dp,2.0_dp]; status=[1,0]; times=y
  call expprep2_summary(tab,x,y,status,times,s1)
  call netwei_summary(tab,x,y,status,times,s2)
  if(maxval(abs(s1%yi-s2%yi))>1.0e-12_dp) error stop 1
  print *, 'test_expprep2: PASS'
end program test_expprep2
