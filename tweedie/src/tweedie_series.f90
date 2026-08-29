! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module tweedie_series_mod
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf, ieee_quiet_nan
use r_compat, only: dp, dgamma, dpois, pgamma
implicit none
private
public :: dtweedie_series, ptweedie_series
public :: dtweedie_logw_smallp, dtweedie_jw_smallp
public :: dtweedie_logv_bigp, dtweedie_kv_bigp
real(dp), parameter :: pi=acos(-1.0_dp)
contains

function dtweedie_logw_smallp(y,phi,power,lo,hi) result(logw)
real(dp),intent(in)::y,phi,power
integer,intent(out),optional::lo,hi
real(dp)::logw,a,a1,r,jmax,jr,cc,wmax,est,m,s,aa
integer::jlo,jhi,j
if(power<=1.0_dp.or.power>=2.0_dp.or.phi<=0.0_dp.or.y<=0.0_dp)then
 logw=ieee_value(0.0_dp,ieee_quiet_nan)
 if(present(lo))lo=0
 if(present(hi))hi=0
 return
end if
a=(2.0_dp-power)/(1.0_dp-power)
a1=1.0_dp-a
r=-a*log(y)+a*log(power-1.0_dp)-a1*log(phi)-log(2.0_dp-power)
jmax=y**(2.0_dp-power)/(phi*(2.0_dp-power))
jr=max(1.0_dp,jmax)
cc=r+a1+a*log(-a)
wmax=a1*jmax
est=wmax
do while(est>wmax-37.0_dp)
 jr=jr+2.0_dp
 est=jr*(cc-a1*log(jr))
end do
jhi=ceiling(jr)
jr=max(1.0_dp,jmax)
est=wmax
do while(est>wmax-37.0_dp.and.jr>=2.0_dp)
 jr=max(1.0_dp,jr-2.0_dp)
 est=jr*(cc-a1*log(jr))
end do
jlo=max(1,floor(jr))
m=-huge(1.0_dp)
do j=jlo,jhi
 aa=r*real(j,dp)-log_gamma(real(j+1,dp))-log_gamma(-a*real(j,dp))
 m=max(m,aa)
end do
s=0.0_dp
do j=jlo,jhi
 aa=r*real(j,dp)-log_gamma(real(j+1,dp))-log_gamma(-a*real(j,dp))
 s=s+exp(aa-m)
end do
logw=log(s)+m
if(present(lo))lo=jlo
if(present(hi))hi=jhi
end function

function dtweedie_jw_smallp(y,phi,power,lo,hi) result(jw)
real(dp),intent(in)::y,phi,power
integer,intent(out),optional::lo,hi
real(dp)::jw,a,a1,r,jmax,jr,cc,wmax,est,m,s,aa
integer::jlo,jhi,j
if(power<=1.0_dp.or.power>=2.0_dp.or.phi<=0.0_dp.or.y<=0.0_dp)then
 jw=ieee_value(0.0_dp,ieee_quiet_nan)
 if(present(lo))lo=0
 if(present(hi))hi=0
 return
end if
a=(2.0_dp-power)/(1.0_dp-power)
a1=1.0_dp-a
r=-a*log(y)+a*log(power-1.0_dp)-a1*log(phi)-log(2.0_dp-power)
jmax=y**(2.0_dp-power)/(phi*(2.0_dp-power))
jr=max(1.0_dp,jmax)
cc=r+a1+a*log(-a)
wmax=a1*jmax
est=wmax
do while(est>wmax-37.0_dp)
 jr=jr+2.0_dp
 est=jr*(cc-a1*log(jr))
end do
jhi=ceiling(jr)
jr=max(1.0_dp,jmax)
est=wmax
do while(est>wmax-37.0_dp.and.jr>=2.0_dp)
 jr=max(1.0_dp,jr-2.0_dp)
 est=jr*(cc-a1*log(jr))
end do
jlo=max(1,floor(jr))
m=-huge(1.0_dp)
do j=jlo,jhi
 aa=r*real(j,dp)-log_gamma(real(j+1,dp))-log_gamma(-a*real(j,dp))+log(real(j,dp))
 m=max(m,aa)
end do
s=0.0_dp
do j=jlo,jhi
 aa=r*real(j,dp)-log_gamma(real(j+1,dp))-log_gamma(-a*real(j,dp))+log(real(j,dp))
 s=s+exp(aa-m)
end do
jw=s*exp(m)
if(present(lo))lo=jlo
if(present(hi))hi=jhi
end function

function dtweedie_logv_bigp(y,phi,power,lo,hi) result(logv)
real(dp),intent(in)::y,phi,power
integer,intent(out),optional::lo,hi
real(dp)::logv,a,a1,r,kmax,kr,cc,vmax,est,m,s,aa,c
integer::klo,khi,k
if(power<=2.0_dp.or.phi<=0.0_dp.or.y<=0.0_dp)then
 logv=ieee_value(0.0_dp,ieee_quiet_nan)
 if(present(lo))lo=0
 if(present(hi))hi=0
 return
end if
a=(2.0_dp-power)/(1.0_dp-power)
a1=1.0_dp-a
r=-a1*log(phi)-log(power-2.0_dp)-a*log(y)+a*log(power-1.0_dp)
kmax=y**(2.0_dp-power)/(phi*(power-2.0_dp))
kr=max(1.0_dp,kmax)
cc=r+a1+a*log(a)
vmax=kmax*a1
est=vmax
do while(est>vmax-37.0_dp)
 kr=kr+2.0_dp
 est=kr*(cc-a1*log(kr))
end do
khi=ceiling(kr)
kr=max(1.0_dp,kmax)
est=vmax
do while(est>vmax-37.0_dp.and.kr>=2.0_dp)
 kr=max(1.0_dp,kr-2.0_dp)
 est=kr*(cc-a1*log(kr))
end do
klo=max(1,floor(kr))
m=-huge(1.0_dp)
do k=klo,khi
 aa=r*real(k,dp)+log_gamma(1.0_dp+a*real(k,dp))-log_gamma(real(k+1,dp))
 m=max(m,aa)
end do
s=0.0_dp
do k=klo,khi
 aa=r*real(k,dp)+log_gamma(1.0_dp+a*real(k,dp))-log_gamma(real(k+1,dp))
 c=sin(-a*pi*real(k,dp))
 if(mod(k,2)==1)c=-c
 s=s+exp(aa-m)*c
end do
if(s<=0.0_dp)then
logv=ieee_value(0.0_dp,ieee_negative_inf)
else
logv=log(s)+m
end if
if(present(lo))lo=klo
if(present(hi))hi=khi
end function

function dtweedie_kv_bigp(y,phi,power,lo,hi) result(kv)
real(dp),intent(in)::y,phi,power
integer,intent(out),optional::lo,hi
real(dp)::kv,a,a1,r,kmax,kr,cc,vmax,est,m,s,aa,c
integer::klo,khi,k
if(power<=2.0_dp.or.phi<=0.0_dp.or.y<=0.0_dp)then
 kv=ieee_value(0.0_dp,ieee_quiet_nan)
 if(present(lo))lo=0
 if(present(hi))hi=0
 return
end if
a=(2.0_dp-power)/(1.0_dp-power)
a1=1.0_dp-a
r=-a1*log(phi)-log(power-2.0_dp)-a*log(y)+a*log(power-1.0_dp)
kmax=y**(2.0_dp-power)/(phi*(power-2.0_dp))
kr=max(1.0_dp,kmax)
cc=r+a1+a*log(a)
vmax=kmax*a1
est=vmax
do while(est>vmax-37.0_dp)
 kr=kr+2.0_dp
 est=kr*(cc-a1*log(kr))
end do
khi=ceiling(kr)
kr=max(1.0_dp,kmax)
est=vmax
do while(est>vmax-37.0_dp.and.kr>=2.0_dp)
 kr=max(1.0_dp,kr-2.0_dp)
 est=kr*(cc-a1*log(kr))
end do
klo=max(1,floor(kr))
m=-huge(1.0_dp)
do k=klo,khi
 aa=r*real(k,dp)+log_gamma(1.0_dp+a*real(k,dp))-log_gamma(real(k+1,dp))+log(real(k,dp))
 m=max(m,aa)
end do
s=0.0_dp
do k=klo,khi
 aa=r*real(k,dp)+log_gamma(1.0_dp+a*real(k,dp))-log_gamma(real(k+1,dp))+log(real(k,dp))
 c=sin(-a*pi*real(k,dp))
 if(mod(k,2)==1)c=-c
 s=s+exp(aa-m)*c
end do
kv=s*exp(m)
if(present(lo))lo=klo
if(present(hi))hi=khi
end function

function dtweedie_series(y,power,mu,phi,lo,hi) result(d)
real(dp),intent(in)::y,power,mu,phi
integer,intent(out),optional::lo,hi
real(dp)::d,logw,logv,tau,lambda,theta,kappa
integer::l,h
if(phi<=0.0_dp.or.mu<=0.0_dp.or.y<0.0_dp.or.power<1.0_dp)then
 d=ieee_value(0.0_dp,ieee_quiet_nan)
 if(present(lo))lo=0
 if(present(hi))hi=0
 return
end if
if(power==1.0_dp)then
 d=dpois(y/phi,lambda=mu/phi)/phi
 l=0
 h=0
else if(power==2.0_dp)then
 d=dgamma(y,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
 l=0
 h=0
else if(y==0.0_dp)then
 if(power<2.0_dp)then
 lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
 d=exp(-lambda)
 else
 d=0.0_dp
 end if
 l=0
 h=0
else if(power<2.0_dp)then
 logw=dtweedie_logw_smallp(y,phi,power,l,h)
 tau=phi*(power-1.0_dp)*mu**(power-1.0_dp)
 lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
 d=exp(-y/tau-lambda-log(y)+logw)
else
 logv=dtweedie_logv_bigp(y,phi,power,l,h)
 theta=mu**(1.0_dp-power)/(1.0_dp-power)
 kappa=mu**(2.0_dp-power)/(2.0_dp-power)
 d=exp((y*theta-kappa)/phi-log(pi*y)+logv)
end if
if(present(lo))lo=l
if(present(hi))hi=h
end function

function ptweedie_series(q,power,mu,phi,iterations) result(cdf)
real(dp),intent(in)::q,power,mu,phi
integer,intent(out),optional::iterations
real(dp)::cdf,lambda,tau,alpha,logfmax,est,nr,term
integer::lo_n,hi_n,n
if(phi<=0.0_dp.or.mu<=0.0_dp.or.power<=1.0_dp.or.power>=2.0_dp)then
 cdf=ieee_value(0.0_dp,ieee_quiet_nan)
 if(present(iterations))iterations=0
 return
end if
if(q<0.0_dp)then
cdf=0.0_dp
if(present(iterations))iterations=0
return
end if
lambda=mu**(2.0_dp-power)/(phi*(2.0_dp-power))
tau=phi*(power-1.0_dp)*mu**(power-1.0_dp)
alpha=(2.0_dp-power)/(1.0_dp-power)
logfmax=-0.5_dp*log(lambda)
est=logfmax
nr=max(lambda,1.0_dp)
do while(est>logfmax-39.0_dp.and.nr>1.0_dp)
 nr=max(1.0_dp,nr-2.0_dp)
 est=-lambda+nr*(log(lambda)-log(nr)+1.0_dp)-0.5_dp*log(nr)
end do
lo_n=max(1,floor(nr))
nr=max(lambda,1.0_dp)
est=logfmax
do while(est>logfmax-39.0_dp)
 nr=nr+1.0_dp
 est=-lambda+nr*(log(lambda)-log(nr)+1.0_dp)-0.5_dp*log(nr)
 if(nr>=1.0e6_dp)exit
end do
hi_n=max(lo_n,min(1000000,ceiling(nr)))
cdf=exp(-lambda)
do n=lo_n,hi_n
 term=dpois(real(n,dp),lambda=lambda)*pgamma(q,shape=-alpha*real(n,dp),rate=1.0_dp/tau)
 cdf=cdf+term
end do
cdf=min(1.0_dp,max(0.0_dp,cdf))
if(present(iterations))iterations=hi_n-lo_n+1
end function

end module tweedie_series_mod
