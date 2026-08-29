! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module tweedie_distribution_mod
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
use r_compat, only: dp, dpois, ppois, qpois, dgamma, pgamma, qgamma, rpois, rgamma, runif1
use tweedie_math, only: tweedie_dev, tweedie_lambda, tweedie_convert, tweedie_convert_result, &
 dinvgauss_tweedie, pinvgauss_tweedie, tweedie_integrand_values
use tweedie_series_mod, only: dtweedie_series, ptweedie_series
use tweedie_interpolation_mod, only: dtweedie_interp, interpolation_available
use tweedie_inversion_mod, only: dtweedie_inversion, ptweedie_inversion
implicit none
private
public :: dtweedie, ptweedie, qtweedie, rtweedie, rtweedie_vec_params
public :: dtweedie_vec, ptweedie_vec, qtweedie_vec
public :: dtweedie_saddle, tweedie_dev, tweedie_lambda, tweedie_convert, tweedie_convert_result
public :: tweedie_integrand_values, dtweedie_series, ptweedie_series, dtweedie_inversion, ptweedie_inversion
real(dp),parameter::pi=acos(-1.0_dp)
contains

function dtweedie(y,mu,phi,power) result(d)
real(dp),intent(in)::y,mu,phi,power
real(dp)::d,xi,xix,rho,dev
if(power<1.0_dp.or.mu<=0.0_dp.or.phi<=0.0_dp)then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(power==1.0_dp)then
d=dpois(y/phi,lambda=mu/phi)
return
end if
if(power==2.0_dp)then
d=dgamma(y,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
return
end if
if(power==3.0_dp)then
d=dinvgauss_tweedie(y,mu,phi)
return
end if
if(y<0.0_dp)then
d=0.0_dp
return
end if
if(y==0.0_dp)then
 if(power<2.0_dp)then
 d=exp(-mu**(2.0_dp-power)/(phi*(2.0_dp-power)))
 else
 d=0.0_dp
 end if
 return
end if
xi=phi*y**(power-2.0_dp)
xix=xi/(1.0_dp+xi)
if(power<=1.1_dp.or.power>10.0_dp)then
 d=dtweedie_series(y,power,mu,phi)
else if(interpolation_available(power,xix))then
 rho=dtweedie_interp(power,xix)
 dev=tweedie_dev(y,mu,power)
 d=rho/(y*sqrt(2.0_dp*pi*xi))*exp(-dev/(2.0_dp*phi))
 if(d<0.0_dp)d=0.0_dp
else
 d=dtweedie_series(y,power,mu,phi)
end if
end function

function ptweedie(q,mu,phi,power) result(p)
real(dp),intent(in)::q,mu,phi,power
real(dp)::p
if(power<1.0_dp.or.mu<=0.0_dp.or.phi<=0.0_dp)then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(power==1.0_dp)then
p=ppois(q/phi,lambda=mu/phi)
else if(power==2.0_dp)then
p=pgamma(q,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
else if(power==3.0_dp)then
p=pinvgauss_tweedie(q,mu,phi)
else if(q<0.0_dp)then
p=0.0_dp
else if(q==0.0_dp.and.power<2.0_dp)then
p=exp(-mu**(2.0_dp-power)/(phi*(2.0_dp-power)))
else if(q==0.0_dp)then
p=0.0_dp
else
p=ptweedie_inversion(q,mu,phi,power)
end if
p=min(1.0_dp,max(0.0_dp,p))
end function

function qtweedie(prob,mu,phi,power) result(q)
real(dp),intent(in)::prob,mu,phi,power
real(dp)::q,qp,qg,start,a,b,fa,fb,fm,m
integer::iter
if(prob<0.0_dp.or.prob>1.0_dp.or.power<1.0_dp.or.mu<=0.0_dp.or.phi<=0.0_dp)then
q=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(prob==0.0_dp)then
q=0.0_dp
return
end if
if(prob==1.0_dp)then
q=ieee_value(0.0_dp,ieee_positive_inf)
return
end if
if(power==1.0_dp)then
q=phi*qpois(prob,lambda=mu/phi)
return
end if
if(power==2.0_dp)then
q=qgamma(prob,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
return
end if
qp=qpois(prob,lambda=mu/phi)
qg=qgamma(prob,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
if(power<2.0_dp)then
 if(prob<=dtweedie(0.0_dp,mu,phi,power))then
 q=0.0_dp
 return
 end if
 start=(qg-phi*qp)*power+(2.0_dp*phi*qp-qg)
else
 start=qg
end if
start=max(start,tiny(1.0_dp))
fa=ptweedie(start,mu,phi,power)-prob
if(fa==0.0_dp)then
q=start
return
end if
if(fa>0.0_dp)then
 b=start
 a=start
 do iter=1,200
   a=0.5_dp*a
   fb=ptweedie(a,mu,phi,power)-prob
   if(fb<0.0_dp)exit
 end do
else
 a=start
 b=start
 do iter=1,200
   b=1.5_dp*(b+2.0_dp)
   fb=ptweedie(b,mu,phi,power)-prob
   if(fb>0.0_dp)exit
 end do
end if
if(fa>0.0_dp)then
 ! b is the original point; fa is its positive value, and a has the negative value.
else
 ! a is the original point and b was expanded upward.
end if
fa=ptweedie(a,mu,phi,power)-prob
fb=ptweedie(b,mu,phi,power)-prob
if(fa>0.0_dp.or.fb<0.0_dp)then
q=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
do iter=1,200
 m=0.5_dp*(a+b)
 fm=ptweedie(m,mu,phi,power)-prob
 if(abs(fm)<1.0e-12_dp.or.abs(b-a)<=1.0e-12_dp*max(1.0_dp,abs(m)))exit
 if(fm>0.0_dp)then
 b=m
 else
 a=m
 end if
end do
q=0.5_dp*(a+b)
end function

function rtweedie(n,mu,phi,power) result(x)
integer,intent(in)::n
real(dp),intent(in)::mu,phi,power
real(dp)::x(max(0,n)),lambda,alpha,gam
integer::i
integer,allocatable::nn(:)
real(dp),allocatable::gg(:)
if(n<=0)return
if(power<1.0_dp.or.mu<=0.0_dp.or.phi<=0.0_dp)then
x=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(power==1.0_dp)then
 nn=rpois(n,mu/phi)
 x=phi*real(nn,dp)
else if(power==2.0_dp)then
 gg=rgamma(n,shape=1.0_dp/phi,scale=phi*mu)
 x=gg
else if(power>2.0_dp)then
 do i=1,n
 x(i)=qtweedie(runif1(),mu,phi,power)
 end do
else
 lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
 alpha=(2.0_dp-power)/(1.0_dp-power)
 gam=phi*(power-1.0_dp)*mu**(power-1.0_dp)
 nn=rpois(n,lambda)
 do i=1,n
   if(nn(i)==0)then
   x(i)=0.0_dp
   else
   gg=rgamma(1,shape=-real(nn(i),dp)*alpha,scale=gam)
   x(i)=gg(1)
   end if
 end do
end if
end function

function rtweedie_vec_params(n,mu,phi,power) result(x)
integer,intent(in)::n
real(dp),intent(in)::mu(:),phi(:),power
real(dp)::x(max(0,n)),one(1)
integer::i,im,ip
if(n<=0)return
if(size(mu)==0.or.size(phi)==0)then
 x=ieee_value(0.0_dp,ieee_quiet_nan)
 return
end if
do i=1,n
 im=1+mod(i-1,size(mu))
 ip=1+mod(i-1,size(phi))
 one=rtweedie(1,mu(im),phi(ip),power)
 x(i)=one(1)
end do
end function rtweedie_vec_params

pure elemental function dtweedie_saddle(y,mu,phi,power,eps) result(d)
real(dp),intent(in)::y,mu,phi,power
real(dp),intent(in),optional::eps
real(dp)::d,ye,e,dev,lambda
e=1.0_dp/6.0_dp
if(present(eps))e=eps
if(phi<=0.0_dp.or.mu<=0.0_dp.or.(power>=1.0_dp.and.y<0.0_dp))then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
ye=y
if(power<2.0_dp)ye=y+e
dev=tweedie_dev(y,mu,power)
d=(2.0_dp*pi*phi*ye**power)**(-0.5_dp)*exp(-dev/(2.0_dp*phi))
if(y==0.0_dp)then
 if(power>=1.0_dp.and.power<2.0_dp)then
 lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
 d=exp(-lambda)
 else
 d=0.0_dp
 end if
end if
end function

function dtweedie_vec(y,mu,phi,power) result(d)
real(dp),intent(in)::y(:),mu(:),phi(:),power
real(dp)::d(size(y))
integer::i
if(size(mu)/=size(y).or.size(phi)/=size(y))then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
do i=1,size(y)
d(i)=dtweedie(y(i),mu(i),phi(i),power)
end do
end function
function ptweedie_vec(q,mu,phi,power) result(p)
real(dp),intent(in)::q(:),mu(:),phi(:),power
real(dp)::p(size(q))
integer::i
if(size(mu)/=size(q).or.size(phi)/=size(q))then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
do i=1,size(q)
p(i)=ptweedie(q(i),mu(i),phi(i),power)
end do
end function
function qtweedie_vec(prob,mu,phi,power) result(q)
real(dp),intent(in)::prob(:),mu(:),phi(:),power
real(dp)::q(size(prob))
integer::i
if(size(mu)/=size(prob).or.size(phi)/=size(prob))then
q=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
do i=1,size(prob)
q(i)=qtweedie(prob(i),mu(i),phi(i),power)
end do
end function

end module tweedie_distribution_mod
