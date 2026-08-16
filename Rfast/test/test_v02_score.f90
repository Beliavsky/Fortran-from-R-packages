program test_v02_score
   use rfast
   implicit none
   real(dp) :: xr(8,3),yb(8),yg(8),yi(8)
   integer :: yc(8),ygm(8),ynb(8),yztp(8),ym(8)
   type(score_result) :: r
   integer :: i,fail
   fail=0
   do i=1,8
      xr(i,1)=real(i,dp);xr(i,2)=(-1.0_dp)**i;xr(i,3)=sin(real(i,dp))
   end do
   yb=[0.12_dp,0.20_dp,0.31_dp,0.42_dp,0.55_dp,0.63_dp,0.74_dp,0.82_dp]
   yg=[0.4_dp,0.8_dp,1.2_dp,1.8_dp,2.3_dp,2.9_dp,3.5_dp,4.2_dp]
   yi=[0.8_dp,1.1_dp,1.4_dp,1.9_dp,2.2_dp,2.8_dp,3.1_dp,3.7_dp]
   yc=[0,1,0,2,1,3,2,4];ygm=[0,1,1,2,2,3,3,4];ynb=[0,0,1,1,2,3,4,6];yztp=[1,1,2,2,3,3,4,5]
   ym=[1,2,3,1,2,3,2,3]
   r=score_betaregs(yb,xr);call check(r,'beta')
   r=score_expregs(yg,xr);call check(r,'exp')
   r=score_gammaregs(yg,xr);call check(r,'gamma')
   r=score_geomregs(ygm,xr);call check(r,'geom')
   r=score_glms(real(yc,dp),xr,.false.);call check(r,'glm')
   r=score_invgaussregs(yi,xr);call check(r,'invgauss')
   r=score_multinomregs(ym,xr);call check(r,'multinom')
   r=score_negbinregs(ynb,xr);call check(r,'negbin')
   r=score_weibregs(yg,xr);call check(r,'weib')
   r=score_ztpregs(yztp,xr);call check(r,'ztp')
   if(fail==0)then;print '(a)','test_v02_score: PASS';else;print '(a,i0)','test_v02_score: FAIL ',fail;error stop 1;end if
contains
   subroutine check(a,name)
      type(score_result),intent(in)::a;character(*),intent(in)::name
      if(.not.allocated(a%statistic).or..not.allocated(a%pvalue))then;print *,name,' allocation';fail=fail+1;return;end if
      if(any(a%statistic<0.0_dp).or.any(a%pvalue<0.0_dp).or.any(a%pvalue>1.0_dp))then;print *,name,' range';fail=fail+1;end if
   end subroutine check
end program test_v02_score
