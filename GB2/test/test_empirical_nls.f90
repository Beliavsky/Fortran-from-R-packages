program test_empirical_nls
  use gb2, only : dp,main_emp,robust_weights,main2_gb2,qgb2,nlsfit_gb2,optimization_result
  implicit none
  real(dp)::x(6),w(6),v(6),corr(6),aw(6),target(6),ei4(4),par0(4),med
  type(optimization_result)::fit
  integer::fails
  x=[1._dp,2._dp,3._dp,4._dp,5._dp,6._dp]
  w=1._dp
  fails=0
  call main_emp(x,w,v)
  if(abs(v(1)-3._dp)>1e-12_dp .or. abs(v(2)-3.5_dp)>1e-12_dp) fails=fails+1
  call robust_weights(x,w,corr=corr,adjusted=aw)
  if(any(corr<=0._dp) .or. any(corr>1._dp) .or. &
    maxval(abs(aw-corr))>1e-14_dp) fails=fails+1
  par0=[2._dp,3.5_dp,2._dp,3._dp]
  call main2_gb2(.6_dp,par0(1),1._dp,par0(1)*par0(3),par0(1)*par0(4),target)
  ei4=target(3:6)
  med=qgb2(.5_dp,par0(1),par0(2),par0(3),par0(4))
  call nlsfit_gb2(med,ei4,par0,fit)
  if(.not.fit%converged .or. maxval(abs(fit%par-par0))>2e-6_dp) then
  print *,'nls fit',fit%par
  fails=fails+1
  end if
  if(fails>0) error stop 1
  print '(a)','test_empirical_nls: PASS'
end program
