program test_v02_repeated
   use rfast
   implicit none
   real(dp)::y(8,3),yb(8,3,2),xv(3),vc(8,2),mom(2,5),rr(2),rma(2,2),rml(2,2)
   integer::id(8),i,fail,idu(11)
   real(dp)::xu(11,1),refinfo(3),refran(4)
   type(variance_components_result)::vm,vu
   fail=0
   do i=1,8
      y(i,:)=[real(i,dp)*0.1_dp,2.0_dp+real(i,dp)*0.1_dp,4.0_dp+real(i,dp)*0.1_dp]
   end do
   rr=rm_anova(y);if(rr(1)<100.0_dp.or.rr(2)>1e-6_dp)fail=fail+1
   yb(:,:,1)=y;yb(:,:,2)=0.5_dp*y+0.2_dp;rma=rm_anovas(yb)
   if(rma(1,1)<=0.0_dp)fail=fail+1
   xv=[-1.0_dp,0.0_dp,1.0_dp];rml=rm_lines(yb,xv)
   if(any(rml(:,2)<0.0_dp))fail=fail+1
   id=[1,1,2,2,3,3,4,4]
   do i=1,8;vc(i,1)=real(id(i),dp)+0.05_dp*(-1.0_dp)**i;vc(i,2)=2.0_dp*real(id(i),dp)+0.1_dp*sin(real(i,dp));end do
   mom=varcomps_mom(vc,id);if(any(mom(:,1)<0.0_dp).or.any(mom(:,2)<=0.0_dp))fail=fail+1
   vm=varcomps_mle_balanced(vc,id,.true.)
   if(vm%status/=0.or.any(vm%info(:,1)<0.0_dp).or..not.allocated(vm%ranef))fail=fail+1
   idu=[1,1,1,2,2,3,3,3,3,4,4]
   xu(:,1)=[0.9_dp,1.1_dp,1.0_dp,2.0_dp,2.2_dp,3.0_dp,3.2_dp,2.9_dp,3.1_dp,4.0_dp,4.1_dp]
   vu=colvarcomps_mle(xu,idu,.true.)
   refinfo=[1.177074772627946_dp,0.0235414954543306_dp,-3.239895251603404_dp]
   refran=[-1.539079614155474_dp,-0.444891232591037_dp,0.498169010032368_dp,1.485801836714455_dp]
   if(vu%status/=0.or.maxval(abs(vu%info(1,:)-refinfo))>1e-10_dp)fail=fail+1
   if(.not.allocated(vu%ranef))then
      fail=fail+1
   else if(maxval(abs(vu%ranef(:,1)-refran))>1e-10_dp)then
      fail=fail+1
   end if
   if(fail==0)then;print '(a)','test_v02_repeated: PASS';else;print '(a,i0)','test_v02_repeated: FAIL ',fail;error stop 1;end if
end program test_v02_repeated
