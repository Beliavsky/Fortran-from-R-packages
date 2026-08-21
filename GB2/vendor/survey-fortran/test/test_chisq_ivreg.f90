program test_chisq_ivreg
  use survey
  implicit none
  type(survey_design_t)::d,d2
  type(chisq_result_t)::cw,cr
  type(ivreg_result_t)::iv
  integer::row(12),col(12),cl(12,1),cl2(8,1),i,fails
  real(dp)::w(12),w2(8),x(8,2),z(8,3),y(8),z1(8),z2(8),u(8),v(8)
  fails=0;row=[1,1,1,1,1,2,2,2,2,2,2,2];col=[1,1,1,2,2,1,1,2,2,2,2,2];w=1
  do i=1,12;cl(i,1)=i;end do;call make_design(w,cl,d)
  call svy_chisq_wald(row,2,col,2,d,cw);call near(cw%statistic,1.159407665505228_dp,2e-10_dp,'Wald F',fails)
  call svy_chisq_rao_scott(row,2,col,2,d,cr);call near(cr%statistic,1.086530612244898_dp,2e-10_dp,'Rao-Scott F',fails);call near(cr%p_value,0.3196089904931203_dp,2e-9_dp,'Rao-Scott p',fails)
  z1=[-2._dp,-1._dp,0._dp,1._dp,2._dp,3._dp,4._dp,5._dp];z2=[1._dp,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,0._dp]
  u=[.2_dp,-.1_dp,.3_dp,-.2_dp,.1_dp,-.3_dp,.2_dp,-.1_dp];v=[.1_dp,.2_dp,-.2_dp,.1_dp,-.1_dp,.3_dp,-.2_dp,.1_dp]
  w2=[1._dp,2._dp,1._dp,2._dp,1._dp,2._dp,1._dp,2._dp]
  do i=1,8;cl2(i,1)=i;x(i,1)=1;x(i,2)=.8_dp*z1(i)+.4_dp*z2(i)+v(i);z(i,:)=[1._dp,z1(i),z2(i)];y(i)=1.2_dp+2.5_dp*x(i,2)+u(i);end do
  call make_design(w2,cl2,d2);call svy_ivreg(x,z,y,d2,iv)
  call near(iv%coef(1),1.1902800297997285_dp,2e-10_dp,'IV intercept',fails);call near(iv%coef(2),2.474012884000175_dp,2e-10_dp,'IV slope',fails)
  if(fails>0)error stop 1;print '(a)','test_chisq_ivreg: PASS'
contains
  subroutine near(a,b,tol,name,f);real(dp),intent(in)::a,b,tol;character(*),intent(in)::name;integer,intent(inout)::f
    if(abs(a-b)>tol*(1+abs(b)))then;print '(a,2es24.15)',trim(name)//' FAIL ',a,b;f=f+1;end if
  end subroutine
end program
