program test_fit
  use flexsurv
  implicit none
  type(flexsurv_data)::dat
  type(flexsurv_spec)::sp
  type(flexsurv_result)::fit
  real(dp)::t(8),rate,ll,manual
  integer::st(8),fails,i
  fails=0;t=[0.2_dp,0.5_dp,0.9_dp,1.1_dp,1.8_dp,2.0_dp,2.5_dp,3.0_dp];st=1
  call prepare_survival_data(dat,t,st)
  call initialize_spec(sp,dist_exponential,size(t),[0.8_dp])
  fit=fit_flexsurvreg(dat,sp)
  rate=real(size(t),dp)/sum(t)
  if(.not.fit%converged)then;print *,'fit did not converge';fails=fails+1;end if
  if(abs(fit%base(1)-rate)>2e-5_dp)then;print *,'rate ',fit%base(1),rate;fails=fails+1;end if
  manual=real(size(t),dp)*log(rate)-rate*sum(t)
  if(abs(fit%loglik-manual)>1e-6_dp)then;print *,'ll ',fit%loglik,manual;fails=fails+1;end if
  ! Regression on rate: two rows with log-rate shift.
  call initialize_spec(sp,dist_exponential,2,[0.5_dp])
  deallocate(sp%reg(1)%x);allocate(sp%reg(1)%x(2,1));sp%reg(1)%x(:,1)=[0.0_dp,1.0_dp]
  if(abs(predict_survival(sp,[log(0.5_dp),log(2.0_dp)],1,1.0_dp)-exp(-0.5_dp))>1e-12_dp)fails=fails+1
  if(abs(predict_survival(sp,[log(0.5_dp),log(2.0_dp)],2,1.0_dp)-exp(-1.0_dp))>1e-12_dp)fails=fails+1
  ! Right censoring / left truncation likelihood identity.
  t=[0.2_dp,0.5_dp,0.9_dp,1.1_dp,1.8_dp,2.0_dp,2.5_dp,3.0_dp];st=[1,0,1,0,1,1,0,1]
  call prepare_survival_data(dat,t,st,start=[(0.05_dp*real(i-1,dp),i=1,8)])
  call initialize_spec(sp,dist_exponential,8,[0.7_dp])
  ll=flexsurv_loglik(dat,sp,[log(0.7_dp)])
  manual=0.0_dp
  do i=1,8
    if(st(i)==1)manual=manual+log(0.7_dp)-0.7_dp*t(i)
    if(st(i)==0)manual=manual-0.7_dp*t(i)
    manual=manual+0.7_dp*dat%start(i)
  end do
  if(abs(ll-manual)>1e-10_dp)then;print *,'cens ll ',ll,manual;fails=fails+1;end if
  if(fails>0)error stop 1
  print *,'test_fit: PASS'
end program test_fit
