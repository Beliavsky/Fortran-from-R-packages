program test_additional
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use statmod_utils
use statmod_models
use statmod_qres
use r_compat, only: dp, set_seed_int
implicit none
real(dp)::a(2,2),mv(2,2),vm(2,2),pval(1),sage(1),sc
real(dp)::yy(4),xx(4,3),ysage(1)
real(dp)::ymat(2,4),lib(4),disp
real(dp),allocatable::coef(:,:),fit(:,:),qr(:)
integer,allocatable::ord(:)
type(reml_fit_result_t)::rg
real(dp)::xg(6,2),zg(6,1),yg(6)
integer::fail
fail=0
a=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2])
mv=matvec(a,[10.0_dp,100.0_dp])
vm=vecmat([10.0_dp,100.0_dp],a)
if(maxval(abs(mv-reshape([10.0_dp,20.0_dp,300.0_dp,400.0_dp],[2,2])))>1e-14_dp)then
print *,'FAIL matvec'
fail=fail+1
end if
if(maxval(abs(vm-reshape([10.0_dp,200.0_dp,30.0_dp,400.0_dp],[2,2])))>1e-14_dp)then
print *,'FAIL vecmat'
fail=fail+1
end if

yy=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
xx(:,1)=[4.0_dp,1.0_dp,3.0_dp,2.0_dp]
xx(:,2)=yy
xx(:,3)=[1.0_dp,-1.0_dp,1.0_dp,-1.0_dp]
ord=forward_select(yy,xx,intercept=.false.,nvar=1)
if(size(ord)/=1.or.ord(1)/=2)then
print *,'FAIL forward'
fail=fail+1
end if
sc=mscale([-2.0_dp,-1.0_dp,-0.5_dp,0.5_dp,1.0_dp,2.0_dp])
if(.not.ieee_is_finite(sc).or.sc<=0)then
print *,'FAIL mscale'
fail=fail+1
end if
pval=permp([0],10,2,2,total_nperm=3_8,method_exact=.true.)
call close(pval(1),((2.0_dp/3.0_dp)**10+(1.0_dp/3.0_dp)**10)/3.0_dp,1e-13_dp,'permp exact')
sage=sage_test([1],[3],4,4)
call close(sage(1),0.625_dp,1e-14_dp,'sage equal libraries')

ymat=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,2.0_dp,1.0_dp,4.0_dp,3.0_dp],[2,4])
lib=100.0_dp
call fit_nbp(ymat,[1,1,2,2],lib,coef,fit,disp)
call close(fit(1,1),fit(1,2),1e-12_dp,'fitNBP equal group fit 1')
call close(fit(1,3),fit(1,4),1e-12_dp,'fitNBP equal group fit 2')
if(.not.ieee_is_finite(disp).or.disp<0.0_dp)then
print *,'FAIL fitNBP dispersion'
fail=fail+1
end if

xg(:,1)=1.0_dp
xg(:,2)=[-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,1.5_dp]
zg=1.0_dp
yg=[0.8_dp,1.0_dp,1.2_dp,1.5_dp,1.9_dp,2.4_dp]
rg=remlscoregamma(yg,xg,zg,maxit=12)
if(.not.allocated(rg%beta).or.any(.not.ieee_is_finite(rg%beta)).or.any(rg%phi<=0))then
 print *,'FAIL remlscoregamma'
 fail=fail+1
end if
call set_seed_int(999)
call qres_tweedie([0.0_dp,1.0_dp],[1.0_dp,1.5_dp],1.5_dp,df=1,dispersion=0.8_dp,resid=qr)
if(any(.not.ieee_is_finite(qr)))then
print *,'FAIL qres tweedie'
fail=fail+1
end if
if(fail>0)error stop 'test_additional failed'
print '(a)','test_additional: PASS'
contains
subroutine close(x,y,tol,name)
real(dp),intent(in)::x,y,tol
character(len=*),intent(in)::name
if(.not.ieee_is_finite(x).or.abs(x-y)>tol*max(1.0_dp,abs(y)))then
 print '(a,2es24.15)','FAIL '//trim(name)//': ',x,y
 fail=fail+1
end if
end subroutine
end program
