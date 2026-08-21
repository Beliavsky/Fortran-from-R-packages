program test_multivariate_pps
  use survey
  implicit none
  type(survey_design_t)::d
  real(dp)::items(6,3),w(6),prob(3),dc(3,3),xx(3,1),tot(1),vv(1,1),alpha
  integer::cl(6,1),i,fails
  fails=0;items=reshape([1._dp,2._dp,3._dp,4._dp,5._dp,6._dp, 2._dp,2._dp,4._dp,4._dp,6._dp,6._dp, 1._dp,2._dp,3._dp,5._dp,5._dp,7._dp],[6,3]);w=1
  do i=1,6;cl(i,1)=i;end do;call make_design(w,cl,d);alpha=svy_cralpha(items,d);call near(alpha,0.9770916334661353_dp,1e-11_dp,'alpha',fails)
  prob=[.5_dp,.6_dp,.7_dp];xx(:,1)=[10._dp,20._dp,30._dp];call poisson_dcheck(prob,dc);call pps_total(xx,prob,dc,tot,vv)
  call near(tot(1),96.1904761904762_dp,1e-11_dp,'PPS total',fails);call near(vv(1,1),1195.46485260771_dp,1e-10_dp,'PPS variance',fails)
  if(fails>0)error stop 1;print '(a)','test_multivariate_pps: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
