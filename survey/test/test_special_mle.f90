module test_mle_callbacks
  use survey, only : dp
  implicit none
  real(dp), save :: yobs(6)
contains
  function llobs(theta,i) result(v)
    real(dp),intent(in)::theta(:);integer,intent(in)::i;real(dp)::v,s
    s=exp(theta(2));v=-log(s)-0.5_dp*((yobs(i)-theta(1))/s)**2
  end function
  subroutine scobs(theta,i,scr)
    real(dp),intent(in)::theta(:);integer,intent(in)::i;real(dp),intent(out)::scr(:);real(dp)::s,r
    s=exp(theta(2));r=yobs(i)-theta(1);scr(1)=r/(s*s);scr(2)=-1+r*r/(s*s)
  end subroutine
end module test_mle_callbacks

program test_special_mle
  use survey
  use test_mle_callbacks
  implicit none
  type(survey_design_t)::d
  type(mle_result_t)::fit
  real(dp)::w(6),start(2)
  integer::cl(6,1),i,fails
  fails=0;call near(chisq_survival(2.5_dp,3._dp),0.4752910833430205_dp,2e-10_dp,'chi-square sf',fails)
  yobs=[1._dp,1.5_dp,2._dp,2.5_dp,3._dp,4._dp];w=[1._dp,2._dp,1._dp,2._dp,1._dp,1._dp];do i=1,6;cl(i,1)=i;end do;call make_design(w,cl,d)
  start=[2._dp,0._dp];call svy_mle(llobs,d,start,fit,score=scobs,maxfun=3000)
  call near(fit%par(1),2.25_dp,3e-5_dp,'MLE mu',fails);call near(fit%par(2),-0.1038196823891223_dp,5e-5_dp,'MLE log sigma',fails)
  if(.not.fit%converged)then;print *,'MLE optimizer status ',fit%status;fails=fails+1;end if
  if(fails>0)error stop 1;print '(a)','test_special_mle: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
