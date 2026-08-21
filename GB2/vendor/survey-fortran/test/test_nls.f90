module test_nls_callbacks
  use survey, only : dp
  implicit none
  real(dp),save::xobs(7)
contains
  function expo(theta,i) result(mu)
    real(dp),intent(in)::theta(:);integer,intent(in)::i;real(dp)::mu
    mu=theta(1)*exp(theta(2)*xobs(i))
  end function
end module
program test_nls
  use survey
  use test_nls_callbacks
  implicit none
  type(survey_design_t)::d;type(nls_result_t)::fit
  real(dp)::y(7),w(7),start(2);integer::cl(7,1),i,fails
  fails=0;xobs=[0._dp,.5_dp,1._dp,1.5_dp,2._dp,2.5_dp,3._dp];y=[1.9_dp,2.6_dp,3.5_dp,4.8_dp,6.4_dp,8.9_dp,11.8_dp];w=[1._dp,2._dp,1._dp,2._dp,1._dp,2._dp,1._dp]
  do i=1,7;cl(i,1)=i;end do;call make_design(w,cl,d);start=[2._dp,.6_dp];call svy_nls(expo,y,d,start,fit)
  call near(fit%coef(1),1.9272539567905491_dp,2e-6_dp,'nls a',fails);call near(fit%coef(2),0.6069451730398547_dp,2e-6_dp,'nls b',fails)
  if(.not.fit%converged)then;print *,'NLS did not converge';fails=fails+1;end if
  if(fails>0)error stop 1;print '(a)','test_nls: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
