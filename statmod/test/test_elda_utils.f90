program test_elda_utils
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use statmod_elda
use statmod_glm_misc
use statmod_utils
use r_compat, only: dp, set_seed_int
implicit none
type(elda_one_group_result_t)::e
real(dp)::response(2),dose(2),tested(2),score(1),r(4),w(4),xe(4,1),xn(4,1)
real(dp)::p(3),pa(3),expected,stat,pv
logical::h(3)
integer::fail
fail=0
response=[1.0_dp,2.0_dp]
dose=1.0_dp
tested=10.0_dp
e=elda_one_group(response,dose,tested)
expected=-log(0.85_dp)
call close(e%lambda,expected,2e-10_dp,'ELDA equal-dose MLE')

r=[-1.0_dp,0.0_dp,1.0_dp,0.0_dp]
w=1.0_dp
xe=1.0_dp
xn(:,1)=[-1.0_dp,0.0_dp,1.0_dp,2.0_dp]
score=glm_scoretest(r,w,xe,xn,dispersion=1.0_dp)
! x_new residualized around mean 0.5: [-1.5,-0.5,0.5,1.5]; numerator=2, ss=5.
call close(score(1),2.0_dp/sqrt(5.0_dp),2e-12_dp,'GLM score')

p=[0.01_dp,0.04_dp,0.20_dp]
pa=p_adjust_holm(p)
call close(pa(1),0.03_dp,1e-14_dp,'Holm 1')
call close(pa(2),0.08_dp,1e-14_dp,'Holm 2')
call close(pa(3),0.20_dp,1e-14_dp,'Holm 3')
h=hommel_test(p)
if(.not.all(h .eqv. [(.false.),(.true.),(.true.)]))then
   ! statmod::hommel.test returns a logical decision vector; this checks its literal algorithm.
   print *,'FAIL Hommel'
   fail=fail+1
end if
call set_seed_int(123)
call compare_two_growth_curves([1,1,2,2],reshape([1.0_dp,1.1_dp,2.0_dp,2.1_dp, &
   2.0_dp,2.1_dp,3.0_dp,3.1_dp],[4,2]),nsim=25,stat=stat,p_value=pv)
if(.not.ieee_is_finite(stat).or.pv<=0.0_dp.or.pv>1.0_dp)then
 print *,'FAIL growth curves'
 fail=fail+1
end if
if(fail>0)error stop 'test_elda_utils failed'
print '(a)','test_elda_utils: PASS'
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
