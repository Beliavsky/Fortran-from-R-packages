program test_models
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use statmod_models
use r_compat, only: dp
implicit none
real(dp)::x1(1,1),y1(1),x5(5,2),y5(5),z7(7,4),x7(7,1),y7(7)
type(glm_fit_result_t)::gg,nb
type(mixed_fit_result_t)::mm
type(reml_fit_result_t)::rr
real(dp)::xr(8,2),zr(8,2),yr(8)
integer::i,fail
fail=0
x1=1
y1=1
gg=glmgam_fit(x1,y1)
call close(gg%coefficients(1),1.0_dp,1e-12_dp,'glmgam coef')
call close(gg%deviance,0.0_dp,1e-12_dp,'glmgam dev')

x5(:,1)=1
x5(:,2)=[1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
y5=[0.0_dp,0.0_dp,6.0_dp,2.0_dp,9.0_dp]
nb=glmnb_fit(x5,y5,[0.1_dp])
call close(nb%coefficients(1),1.734601_dp,2e-5_dp,'glmnb intercept')
call close(nb%fitted(3),5.666667_dp,2e-6_dp,'glmnb fitted')
call close(nb%deviance,3.242349_dp,3e-6_dp,'glmnb deviance')

y7=[-1.0_dp,1.0_dp,-2.0_dp,2.0_dp,0.5_dp,1.7_dp,-0.1_dp]
x7=1
z7=0
z7(1:2,1)=1
z7(3:4,2)=1
z7(5:6,3)=1
z7(7,4)=1
mm=mixed_model2_fit(y7,x7,z7)
call close(mm%varcomp(1),2.923462_dp,4e-6_dp,'mixed residual')
call close(mm%varcomp(2),-1.098564_dp,4e-6_dp,'mixed block')
call close(mm%coefficients(1),0.3376358_dp,4e-6_dp,'mixed coef')

xr(:,1)=1
xr(:,2)=[-1.5_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp]
zr(:,1)=1
zr(:,2)=xr(:,2)
yr=[-1.0_dp,-0.3_dp,-0.1_dp,0.4_dp,0.5_dp,1.2_dp,1.4_dp,2.2_dp]
rr=remlscore(yr,xr,zr,maxit=20)
if(.not.allocated(rr%beta).or.any(.not.ieee_is_finite(rr%beta)).or.any(rr%phi<=0))then
 print *,'FAIL remlscore finite'
 fail=fail+1
end if
if(fail>0)error stop 'test_models failed'
print '(a)','test_models: PASS'
contains
subroutine close(a,b,tol,name)
real(dp),intent(in)::a,b,tol
character(len=*),intent(in)::name
if(.not.ieee_is_finite(a).or.abs(a-b)>tol*max(1.0_dp,abs(b)))then
 print '(a,2es24.15)','FAIL '//trim(name)//': ',a,b
 fail=fail+1
end if
end subroutine
end program
