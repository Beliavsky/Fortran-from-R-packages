program test_rtrunc
  use flexsurv
  implicit none
  type(survrtrunc_result)::np
  real(dp)::t(5),rt(5),ti(5),ll,manual,shape,rate,theta,y
  integer::fails,i
  fails=0;t=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp];rt=20.0_dp;ti=0.0_dp
  np=survrtrunc_fit(t,rt,20.0_dp)
  if(size(np%time)/=5.or.abs(np%surv(5))>1e-14_dp)fails=fails+1
  shape=2.0_dp;rate=3.0_dp;theta=0.2_dp
  ll=flexsurvrtrunc_loglik(t,ti,rt,20.0_dp,dist_gamma,[shape,rate],theta,rtrunc_final)
  manual=0.0_dp
  do i=1,5
    y=ti(i)+t(i)
    manual=manual+dgamma_fs(t(i),shape,rate+theta,.true.) &
      -log(pgamma_fs(y,shape,rate+theta))
  end do
  if(abs(ll-manual)>2e-10_dp)then;print *,'rtrunc ll ',ll,manual;fails=fails+1;end if
  if(fails>0)error stop 1
  print *,'test_rtrunc: PASS'
end program test_rtrunc
