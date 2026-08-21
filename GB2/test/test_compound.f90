program test_compound
  use gb2, only : dp,fg_cgb2,dcgb2,pcgb2,component_moments_cgb2,moment_cgb2,qcgb2, &
    vofp_cgb2,pofv_cgb2,scoreu_cgb2,scores_cgb2,hess_cgb2,fit_cgb2,optimization_result
  implicit none
  real(dp)::pl0(3),pl(3),fac(1,3),facl(1,3),m(3),v(2),plback(3),u(6,2),ff(6,3),g(2),gp(2),gm(2),h(2,2),hn(2,2),vp(2),vm(2),eps
  integer::fails,i
  type(optimization_result)::fit
  pl0=[.3_dp,.4_dp,.3_dp]
  pl=[.2_dp,.5_dp,.3_dp]
  fails=0
  call fg_cgb2([3.1_dp],2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,fac)
  if(maxval(abs(fac(1,:)-[.7833974492763923_dp,1.1970240902193738_dp,.9539037637644423_dp]))>3e-12_dp) fails=fails+1
  call fg_cgb2([3.1_dp],2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,facl,'l')
  if(maxval(abs(facl(1,:)-[.3749685029389472_dp,1.5517044881811892_dp,.8894255128194666_dp]))>4e-12_dp) fails=fails+1
  if(abs(dcgb2(3.1_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,pl)-.31006994617102884_dp)>2e-12_dp) fails=fails+1
  if(abs(pcgb2(3.1_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,pl)-.5692172100386271_dp)>3e-9_dp) fails=fails+1
  call component_moments_cgb2(1._dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,m)
  if(maxval(abs(m-[4.257301863874457_dp,3.034773186512535_dp,2.356466761570338_dp]))>5e-12_dp) fails=fails+1
  if(abs(moment_cgb2(1._dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,pl)-3.07578699450226_dp)>5e-12_dp) fails=fails+1
  if(abs(qcgb2(.5_dp,2.3_dp,4.2_dp,1.7_dp,3.4_dp,pl0,pl)-2.8841662769052054_dp)>5e-7_dp) fails=fails+1
  call vofp_cgb2(pl,v)
  call pofv_cgb2(v,plback)
  if(maxval(abs(pl-plback))>1e-14_dp) fails=fails+1
  ff=reshape([2._dp,.4_dp,1.1_dp,.3_dp,1.7_dp,.8_dp, &
    1.4_dp,.6_dp,1.2_dp,.5_dp,2._dp,.9_dp, &
    1.8_dp,.7_dp,1.1_dp,.6_dp,1.4_dp,1.3_dp],[6,3])
  call scoreu_cgb2(ff,pl,u)
  call scores_cgb2(ff,pl,g)
  call hess_cgb2(u,pl,h)
  call vofp_cgb2(pl,v)
  eps=1e-6_dp
  do i=1,2
    vp=v
    vm=v
    vp(i)=vp(i)+eps
    vm(i)=vm(i)-eps
    call pofv_cgb2(vp,plback)
    call scores_cgb2(ff,plback,gp)
    call pofv_cgb2(vm,plback)
    call scores_cgb2(ff,plback,gm)
    hn(:,i)=(gp-gm)/(2*eps)
  end do
  if(maxval(abs(hn-h/6._dp))>3e-7_dp) then
  print *,'mixture Hessian error',maxval(abs(hn-h/6._dp))
  fails=fails+1
  end if
  call fit_cgb2(ff,pl0,fit,maxiter=500,tol=1e-10_dp)
  call scores_cgb2(ff,fit%par,g)
  if(.not.fit%converged .or. abs(sum(fit%par)-1._dp)>1e-12_dp .or. maxval(abs(g))>2e-8_dp) fails=fails+1
  if(fails>0) error stop 1
  print '(a)','test_compound: PASS'
end program
