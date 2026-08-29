program test_core
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use statmod_invgauss
use statmod_gaussquad
use statmod_expected_deviance
use statmod_glm_misc
use r_compat, only: dp
implicit none
type(quad_rule_t)::g
real(dp)::m,v,q,p
integer::fail
fail=0
call check_close(pinvgauss(0.1_dp,mean=2.0_dp,dispersion=0.5_dp),2.057306477e-5_dp,2e-13_dp,'pinvgauss 1')
call check_close(pinvgauss(1.0_dp,mean=3.0_dp,dispersion=0.5_dp),0.2854596328_dp,2e-10_dp,'pinvgauss 2')
call check_close(pinvgauss(3.1_dp,mean=1.0_dp,dispersion=0.5_dp),0.9812161963_dp,2e-10_dp,'pinvgauss 3')
q=qinvgauss(0.001_dp,mean=1.3_dp,dispersion=0.6_dp)
call check_close(q,0.1271035164_dp,5e-9_dp,'qinvgauss .001')
call check_close(qinvgauss(0.5_dp,mean=1.3_dp,dispersion=0.6_dp),0.9446753861_dp,5e-9_dp,'qinvgauss .5')
call check_close(qinvgauss(0.999_dp,mean=1.3_dp,dispersion=0.6_dp),9.2602074131_dp,5e-8_dp,'qinvgauss .999')

g=gauss_quad(5,'legendre')
call check_close(g%nodes(1),-0.9061798459_dp,2e-10_dp,'legendre node')
call check_close(g%weights(3),0.5688888889_dp,2e-10_dp,'legendre weight')
g=gauss_quad(5,'jacobi',5.0_dp,1.1_dp)
call check_close(g%nodes(1),-0.8844049819_dp,3e-10_dp,'jacobi node')
call check_close(g%weights(1),0.40981005618_dp,3e-10_dp,'jacobi weight')
g=gauss_quad_prob(5,'normal')
call check_close(g%nodes(1),-2.856970014_dp,3e-9_dp,'normal prob node')
call check_close(g%weights(3),0.53333333333_dp,3e-10_dp,'normal prob weight')

call expected_deviance_binomial(0.4_dp,2,m,v)
call check_close(m,1.361204081_dp,8e-9_dp,'binom Edev mean')
call check_close(v,1.802700721_dp,8e-9_dp,'binom Edev var')
call expected_deviance(1.0_dp,'Gamma',2.0_dp,m,v)
call check_close(m,1.081451382_dp,8e-9_dp,'gamma Edev mean')
call check_close(v,2.31894507_dp,8e-8_dp,'gamma Edev var')
call expected_deviance(1.0_dp,'negative.binomial',2.0_dp,m,v)
call check_close(m,1.057480184_dp,8e-9_dp,'NB Edev mean')
call check_close(v,0.9740485644_dp,8e-9_dp,'NB Edev var')
call expected_deviance(2.0_dp,'negative.binomial',2.0_dp,m,v)
call check_close(m,1.120623536_dp,8e-9_dp,'NB2 Edev mean')
call check_close(v,1.6273323121_dp,8e-9_dp,'NB2 Edev var')
call expected_deviance(2.0_dp,'poisson',m=m,v=v)
call check_close(m,1.139404056_dp,8e-9_dp,'Poisson Edev mean')
call check_close(v,2.232975219_dp,8e-9_dp,'Poisson Edev var')

call check_close(tweedie_linkinv(log(2.5_dp),0.0_dp),2.5_dp,1e-14_dp,'Tweedie link inverse')
call check_close(tweedie_variance(2.0_dp,1.25_dp),2.0_dp**1.25_dp,1e-14_dp,'Tweedie variance')
if(fail>0)error stop 'test_core failed'
print '(a)','test_core: PASS'
contains
subroutine check_close(a,b,tol,name)
real(dp),intent(in)::a,b,tol
character(len=*),intent(in)::name
if(.not.ieee_is_finite(a).or.abs(a-b)>tol*max(1.0_dp,abs(b)))then
 print '(a,2es24.15)','FAIL '//trim(name)//': ',a,b
 fail=fail+1
end if
end subroutine
end program
