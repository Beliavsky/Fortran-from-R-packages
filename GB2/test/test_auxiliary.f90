program test_auxiliary
  use gb2, only : dp,pkl_cavgb2,loglik_cavgb2,scores_cavgb2,scoreu_cavgb2,scorez_cavgb2, &
    hess_cavgb2,fit_cavgb2,optimization_result
  implicit none
  real(dp)::fac(8,3),z(8,2),lam(2,2),p(8,3),g(2,2),gn(2,2),lp(2,2),lm(2,2),hstep,u(8,2),sc(8,4),h(4,4),hn(4,4),gp(2,2),gm(2,2)
  integer::i,j,fails
  type(optimization_result)::fit
  fac=reshape([ &
    1.8_dp,.5_dp,1.2_dp,.7_dp,1.6_dp,.6_dp,1.4_dp,.8_dp, &
    .5_dp,1.7_dp,.8_dp,1.5_dp,.6_dp,1.8_dp,.7_dp,1.4_dp, &
    1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp],[8,3])
  z(:,1)=1._dp
  z(:,2)=[-1.5_dp,-1._dp,-.5_dp,0._dp,.5_dp,1._dp,1.5_dp,2._dp]
  lam=reshape([.2_dp,-.1_dp,-.3_dp,.25_dp],[2,2])
  call pkl_cavgb2(z,lam,p)
  fails=0
  if(maxval(abs(sum(p,dim=2)-1._dp))>1e-14_dp .or. minval(p)<=0._dp) fails=fails+1
  call scores_cavgb2(fac,z,lam,g)
  do j=1,2
  do i=1,2
    hstep=1e-6_dp
    lp=lam
    lm=lam
    lp(i,j)=lp(i,j)+hstep
    lm(i,j)=lm(i,j)-hstep
    gn(i,j)=(loglik_cavgb2(fac,z,lp)-loglik_cavgb2(fac,z,lm))/(2*hstep)
  end do
  end do
  if(maxval(abs(g-gn))>2e-7_dp) then
  print *,'aux gradient error',maxval(abs(g-gn))
  fails=fails+1
  end if
  call scoreu_cavgb2(fac,z,lam,u,p)
  call scorez_cavgb2(u,z,sc)
  call hess_cavgb2(u,p,z,h)
  hstep=1e-6_dp
  do j=1,2
  do i=1,2
    lp=lam
    lm=lam
    lp(i,j)=lp(i,j)+hstep
    lm(i,j)=lm(i,j)-hstep
    call scores_cavgb2(fac,z,lp,gp)
    call scores_cavgb2(fac,z,lm,gm)
    hn(:,(j-1)*2+i)=reshape((gp-gm)/(2*hstep),[4])
  end do
  end do
  if(maxval(abs(hn-h/8._dp))>3e-7_dp) then
  print *,'aux Hessian error',maxval(abs(hn-h/8._dp))
  fails=fails+1
  end if
  call fit_cavgb2(fac,z,lam,fit,maxiter=500,tol=1e-9_dp)
  if(.not.fit%converged) then
  print *,'aux fit did not converge'
  fails=fails+1
  end if
  if(fails>0) error stop 1
  print '(a)','test_auxiliary: PASS'
end program
