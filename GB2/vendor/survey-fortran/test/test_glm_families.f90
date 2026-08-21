program test_glm_families
  use survey
  implicit none
  type(survey_design_t)::d
  type(glm_result_t)::bfit,pfit
  real(dp)::xx(10,2),yb(10),yp(10),w(10),xv(10)
  integer::cl(10,1),i,fails
  fails=0;xv=[-2._dp,-1.5_dp,-1._dp,-.5_dp,0._dp,.5_dp,1._dp,1.5_dp,2._dp,2.5_dp]
  yb=[0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,1._dp,1._dp,1._dp]
  yp=[0._dp,1._dp,0._dp,2._dp,1._dp,3._dp,2._dp,5._dp,4._dp,7._dp]
  w=[1._dp,2._dp,1._dp,1._dp,2._dp,1._dp,2._dp,1._dp,1._dp,2._dp]
  do i=1,10;xx(i,:)=[1._dp,xv(i)];cl(i,1)=i;end do;call make_design(w,cl,d)
  call svy_glm(xx,yb,d,bfit,FAMILY_BINOMIAL,LINK_LOGIT);call near(bfit%coef(1),-2.7159694114435395_dp,2e-8_dp,'binomial intercept',fails);call near(bfit%coef(2),2.664917856433701_dp,2e-8_dp,'binomial slope',fails)
  call svy_glm(xx,yp,d,pfit,FAMILY_POISSON,LINK_LOG);call near(pfit%coef(1),0.3758554001404402_dp,2e-8_dp,'poisson intercept',fails);call near(pfit%coef(2),0.6165704308961771_dp,2e-8_dp,'poisson slope',fails)
  if(.not.bfit%converged.or..not.pfit%converged)then;print *,'GLM family fit did not converge';fails=fails+1;end if
  if(fails>0)error stop 1;print '(a)','test_glm_families: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
