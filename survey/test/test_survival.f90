program test_survival
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use survey
  implicit none
  type(survey_design_t)::d
  type(survival_curve_t)::km
  type(cox_result_t)::cx
  type(aft_result_t)::af
  real(dp)::time(8),w(8),x(8,1)
  integer::status(8),cl(8,1),i,fails
  fails=0;time=[1._dp,2._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp];status=[1,1,0,1,1,0,1,1];w=1
  x(:,1)=[-1._dp,0._dp,1._dp,-.5_dp,.5_dp,1.5_dp,-1.5_dp,.8_dp];do i=1,8;cl(i,1)=i;end do;call make_design(w,cl,d)
  call svy_km(time,status,d,km,se=.false.);call near(km%survival(1),1._dp,1e-13_dp,'KM start',fails);call near(km%survival(2),0.875_dp,1e-13_dp,'KM t1',fails)
  call svy_coxph(time,status,x,d,cx);if(.not.cx%converged)then;print *,'cox did not converge';fails=fails+1;end if
  if(.not.all(ieee_is_finite(cx%vcov)))then;print *,'cox vcov nonfinite';fails=fails+1;end if
  call svy_survreg(time,status,x,d,'weibull',af);if(.not.af%converged)then;print *,'survreg did not converge';fails=fails+1;end if
  if(af%scale<=0)then;print *,'survreg nonpositive scale';fails=fails+1;end if
  if(fails>0)error stop 1;print '(a)','test_survival: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
