program test_quantiles_glm
  use survey
  implicit none
  type(survey_design_t)::d
  type(glm_result_t)::fit
  real(dp)::xq(5),wq(5),w(6),xx(6,2),y(6)
  integer::cl(6,1),i,fails
  real(dp)::q
  fails=0;xq=[1._dp,2._dp,3._dp,4._dp,5._dp];wq=1._dp
  q=weighted_quantile(xq,wq,0.25_dp,QRULE_HF7);call near(q,2._dp,1e-13_dp,'HF7 q25',fails)
  q=weighted_quantile(xq,wq,0.40_dp,QRULE_MATH);call near(q,2._dp,1e-13_dp,'math q40',fails)
  q=weighted_quantile(xq,wq,0.40_dp,QRULE_SCHOOL);call near(q,2.5_dp,1e-13_dp,'school q40',fails)
  q=student_t_quantile(0.975_dp,10);call near(q,2.2281388519649385_dp,1e-11_dp,'t q975',fails)
  q=student_t_quantile(0.975_dp,0);call near(q,1.959963984540054_dp,2e-13_dp,'t normal fallback',fails)
  w=[1._dp,2._dp,1._dp,2._dp,1._dp,2._dp];do i=1,6;xx(i,1)=1;xx(i,2)=real(i-1,dp);cl(i,1)=i;end do
  y=[1.6_dp,3.4_dp,5.7_dp,7.3_dp,9.55_dp,11.45_dp];call make_design(w,cl,d);call svy_glm(xx,y,d,fit,FAMILY_GAUSSIAN,LINK_IDENTITY)
  call near(fit%coef(1),1.498717948717949_dp,1e-10_dp,'glm intercept',fails);call near(fit%coef(2),1.985897435897436_dp,1e-10_dp,'glm slope',fails)
  if(.not.fit%converged)then;print *,'glm did not converge';fails=fails+1;end if
  if(fails>0)error stop 1;print '(a)','test_quantiles_glm: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
