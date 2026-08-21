program test_fit
  use gb2, only : dp,qgb2,fit_gb2_full,fit_gb2_profile,optimization_result,scores_gb2
  implicit none
  integer,parameter::n=400
  real(dp)::x(n),w(n),prob,g(4)
  integer::i,fails
  type(optimization_result)::f,p
  do i=1,n
  prob=(real(i,dp)-0.5_dp)/real(n,dp)
  x(i)=qgb2(prob,2.3_dp,4.2_dp,1.7_dp,3.4_dp)
  w(i)=1._dp
  end do
  call fit_gb2_full(x,f,w,maxiter=1000,tol=1e-10_dp)
  call fit_gb2_profile(x,p,w,maxiter=1000,tol=1e-10_dp)
  fails=0
  if(.not.f%converged .or. .not.p%converged) fails=fails+1
  if(maxval(abs(f%par-p%par))>3e-7_dp) then
  print *,'full/profile mismatch',f%par,p%par
  fails=fails+1
  end if
  call scores_gb2(x,f%par,g,w)
  if(maxval(abs(g))>2e-8_dp) then
  print *,'score',g
  fails=fails+1
  end if
  if(fails>0) error stop 1
  print '(a)','test_fit: PASS'
end program
