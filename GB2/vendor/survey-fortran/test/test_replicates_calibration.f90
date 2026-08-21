program test_replicates_calibration
  use survey
  implicit none
  type(rep_design_t)::r
  type(svystat_t)::m
  real(dp)::w(4),y(4,1),mm(6,2),wc(6),neww(6),target(2),ach(2),post(6),popt(2)
  integer::psu(4),cat(6),i,fails
  logical::ok
  fails=0;w=[1.0_dp,2.0_dp,1.0_dp,2.0_dp];y(:,1)=[1.0_dp,2.0_dp,4.0_dp,8.0_dp];psu=[1,2,3,4]
  call make_jk1(psu,w,r);m=rep_mean(y,r);call near(m%variance(1,1),3.931875_dp,1e-12_dp,'JK1 mean variance',fails)
  wc=[1.0_dp,2.0_dp,1.0_dp,2.0_dp,1.0_dp,1.0_dp]
  do i=1,6;mm(i,1)=1.0_dp;mm(i,2)=real(i-1,dp);end do
  target=[12.0_dp,34.0_dp];call calibrate_weights(mm,wc,target,neww,CAL_LINEAR,converged=ok)
  ach=matmul(transpose(mm),neww);if(.not.ok)then;print *,'linear calibration did not converge';fails=fails+1;end if
  call near(ach(1),target(1),1e-10_dp,'cal total 1',fails);call near(ach(2),target(2),1e-10_dp,'cal total 2',fails)
  cat=[1,1,1,2,2,2];popt=[7.0_dp,9.0_dp];call poststratify_weights(cat,2,popt,wc,post)
  call near(sum(post,mask=cat==1),7.0_dp,1e-12_dp,'poststratum 1',fails);call near(sum(post,mask=cat==2),9.0_dp,1e-12_dp,'poststratum 2',fails)
  if(fails>0)error stop 1;print '(a)','test_replicates_calibration: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
