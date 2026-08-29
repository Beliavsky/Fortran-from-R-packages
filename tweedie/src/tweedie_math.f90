! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module tweedie_math
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
use r_compat, only: dp, normal_cdf
implicit none
private
public :: tweedie_dev, tweedie_lambda, tweedie_convert_result, tweedie_convert
public :: dinvgauss_tweedie, pinvgauss_tweedie, tweedie_integrand_values
real(dp), parameter :: pi=acos(-1.0_dp)

type :: tweedie_convert_result
 real(dp) :: poisson_lambda
 real(dp) :: gamma_shape
 real(dp) :: gamma_scale
 real(dp) :: p0
 real(dp) :: gamma_mean
 real(dp) :: gamma_phi
end type
contains

pure elemental function tweedie_dev(y,mu,power) result(dev)
real(dp),intent(in)::y,mu,power
real(dp)::dev,d
if(mu<=0.0_dp)then
dev=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(power==1.0_dp)then
 if(y<0.0_dp)then
 dev=ieee_value(0.0_dp,ieee_quiet_nan)
 else if(y==0.0_dp)then
 dev=2.0_dp*mu
 else
 dev=2.0_dp*(y*log(y/mu)-(y-mu))
 end if
else if(power==2.0_dp)then
 if(y<=0.0_dp)then
 dev=ieee_value(0.0_dp,ieee_positive_inf)
 else
 dev=2.0_dp*(log(mu/y)+y/mu-1.0_dp)
 end if
else if(power==0.0_dp)then
 dev=(y-mu)**2
else
 if(y==0.0_dp .and. power<2.0_dp)then
   d=-(y*mu**(1.0_dp-power))/(1.0_dp-power)+mu**(2.0_dp-power)/(2.0_dp-power)
 else if(y<0.0_dp.and.power>=1.0_dp)then
   dev=ieee_value(0.0_dp,ieee_quiet_nan)
   return
 else
   d=y**(2.0_dp-power)/((1.0_dp-power)*(2.0_dp-power)) &
     -y*mu**(1.0_dp-power)/(1.0_dp-power)+mu**(2.0_dp-power)/(2.0_dp-power)
 end if
 dev=2.0_dp*d
end if
end function

pure elemental function tweedie_lambda(mu,phi,power) result(lambda)
real(dp),intent(in)::mu,phi,power
real(dp)::lambda
if(phi<=0.0_dp.or.mu<=0.0_dp.or.power>=2.0_dp)then
 lambda=ieee_value(0.0_dp,ieee_quiet_nan)
else
 lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
end if
end function

pure function tweedie_convert(mu,phi,power) result(out)
real(dp),intent(in)::mu,phi,power
type(tweedie_convert_result)::out
if(power<=1.0_dp.or.power>=2.0_dp.or.phi<=0.0_dp.or.mu<=0.0_dp)then
 out%poisson_lambda=ieee_value(0.0_dp,ieee_quiet_nan)
 out%gamma_shape=out%poisson_lambda
 out%gamma_scale=out%poisson_lambda
 out%p0=out%poisson_lambda
 out%gamma_mean=out%poisson_lambda
 out%gamma_phi=out%poisson_lambda
 return
end if
out%poisson_lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
out%gamma_shape=(2.0_dp-power)/(power-1.0_dp)
out%gamma_scale=phi*(power-1.0_dp)*mu**(power-1.0_dp)
out%p0=exp(-out%poisson_lambda)
out%gamma_phi=(2.0_dp-power)*(power-1.0_dp)*phi**2*mu**(2.0_dp*(power-1.0_dp))
out%gamma_mean=out%gamma_scale/phi
end function

pure elemental function dinvgauss_tweedie(x,mean,dispersion,log_) result(d)
real(dp),intent(in)::x,mean,dispersion
logical,intent(in),optional::log_
real(dp)::d,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
if(mean<=0.0_dp.or.dispersion<=0.0_dp)then
 d=ieee_value(0.0_dp,ieee_quiet_nan)
 return
end if
if(x<=0.0_dp)then
 if(lg)then
 d=-ieee_value(0.0_dp,ieee_positive_inf)
 else
 d=0.0_dp
 end if
 return
end if
ld=-0.5_dp*log(2.0_dp*pi*dispersion*x**3)-(x-mean)**2/(2.0_dp*dispersion*mean**2*x)
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure elemental function pinvgauss_tweedie(q,mean,dispersion) result(p)
real(dp),intent(in)::q,mean,dispersion
real(dp)::p,lambda,z1,z2,term2
if(mean<=0.0_dp.or.dispersion<=0.0_dp)then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(q<=0.0_dp)then
p=0.0_dp
return
end if
lambda=1.0_dp/dispersion
z1=sqrt(lambda/q)*(q/mean-1.0_dp)
z2=-sqrt(lambda/q)*(q/mean+1.0_dp)
if(2.0_dp*lambda/mean<700.0_dp)then
 term2=exp(2.0_dp*lambda/mean)*normal_cdf(z2)
else
 ! Asymptotic log-Phi approximation avoids overflow in exp(2 lambda/mu).
 term2=exp(2.0_dp*lambda/mean-0.5_dp*z2*z2-log(max(-z2,1.0_dp))-0.5_dp*log(2.0_dp*pi))
end if
p=min(1.0_dp,max(0.0_dp,normal_cdf(z1)+term2))
end function

pure subroutine tweedie_integrand_values(y,power,mu,phi,t,pdf,real_k,imag_k,integrand)
real(dp),intent(in)::y,power,mu,phi,t
logical,intent(in)::pdf
real(dp),intent(out)::real_k,imag_k,integrand
real(dp)::front,alpha,omega
front=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
alpha=(2.0_dp-power)/(1.0_dp-power)
omega=atan((1.0_dp-power)*t*phi/mu**(1.0_dp-power))
real_k=front*(cos(omega*alpha)/(cos(omega)**alpha)-1.0_dp)
imag_k=front*sin(omega*alpha)/(cos(omega)**alpha)-t*y
if(pdf)then
 integrand=exp(real_k)*cos(imag_k)
else if(t==0.0_dp)then
 integrand=mu-y
else
 integrand=exp(real_k)*sin(imag_k)/t
end if
end subroutine

end module tweedie_math
