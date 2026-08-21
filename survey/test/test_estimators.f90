program test_estimators
  use survey
  implicit none
  type(survey_design_t) :: d,dfpc
  type(svystat_t) :: t,m
  real(dp) :: w(4),y(4,1),pop(4,1)
  integer :: cl(4,1)
  integer :: fails
  fails=0;w=[1.0_dp,2.0_dp,1.0_dp,2.0_dp];y(:,1)=[1.0_dp,2.0_dp,4.0_dp,8.0_dp];cl(:,1)=[1,2,3,4]
  call make_design(w,cl,d);t=svy_total(y,d);m=svy_mean(y,d)
  call near(t%estimate(1),25.0_dp,1e-12_dp,'total',fails)
  call near(t%variance(1,1),177.0_dp,1e-11_dp,'total variance',fails)
  call near(m%estimate(1),25.0_dp/6.0_dp,1e-12_dp,'mean',fails)
  call near(m%variance(1,1),3.244855967078189_dp,1e-12_dp,'mean variance',fails)
  pop=8.0_dp;call make_design(w,cl,dfpc,pop_size=pop);t=svy_total(y,dfpc)
  call near(t%variance(1,1),88.5_dp,1e-11_dp,'FPC variance',fails)
  if(fails>0)error stop 1
  print '(a)','test_estimators: PASS'
contains
  subroutine near(a,b,tol,name,f)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
