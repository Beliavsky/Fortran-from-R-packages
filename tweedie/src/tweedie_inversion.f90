! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module tweedie_inversion_mod
use, intrinsic :: iso_c_binding, only: c_int, c_double
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use r_compat, only: dp, dpois, dgamma, ppois, pgamma
use tweedie_math, only: tweedie_dev, dinvgauss_tweedie, pinvgauss_tweedie
implicit none
private
public :: dtweedie_inversion, ptweedie_inversion
interface
 subroutine twcomputation_loop(n,p,phi,y,mu,verbose,pdf,funvalue,exitstatus,relerr,int_regions)
  import :: c_int,c_double
  integer(c_int),intent(in)::n,verbose,pdf
  real(c_double),intent(in)::p,phi(n),y(n),mu(n)
  real(c_double),intent(out)::funvalue(n),relerr(n)
  integer(c_int),intent(out)::exitstatus(n),int_regions(n)
 end subroutine
end interface
contains

function dtweedie_inversion(y,mu,phi,power,method,verbose,exitstatus,relerr,regions,igexact) result(d)
real(dp),intent(in)::y,mu,phi,power
integer,intent(in),optional::method
logical,intent(in),optional::verbose,igexact
integer,intent(out),optional::exitstatus,regions
real(dp),intent(out),optional::relerr
real(dp)::d,dev,theta,kappa,m1,m2,m3,scale
real(c_double)::ya(1),mua(1),phia(1),fv(1),re(1),pp
integer(c_int)::es(1),its(1),n,vb,pdf,im
logical::verb,ig
if(phi<=0.0_dp.or.mu<=0.0_dp.or.power<1.0_dp)then
d=ieee_value(0.0_dp,ieee_quiet_nan)
call set_details(-1,0.0_dp,0)
return
end if
verb=.false.
if(present(verbose))verb=verbose
ig=.true.
if(present(igexact))ig=igexact
if(power==1.0_dp)then
d=dpois(y/phi,lambda=mu/phi)
call set_details(0,0.0_dp,0)
return
end if
if(power==2.0_dp)then
d=dgamma(y,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
call set_details(0,0.0_dp,0)
return
end if
if(power==3.0_dp.and.ig)then
d=dinvgauss_tweedie(y,mu,phi)
call set_details(0,0.0_dp,0)
return
end if
if(y<0.0_dp)then
d=0.0_dp
call set_details(0,0.0_dp,0)
return
end if
if(y==0.0_dp)then
 if(power<2.0_dp)then
 d=exp(-mu**(2.0_dp-power)/(phi*(2.0_dp-power)))
 else
 d=0.0_dp
 end if
 call set_details(0,0.0_dp,0)
 return
end if
! Dunn-Smyth (2008) scaling choice. If method is absent or invalid, choose the smallest multiplier.
theta=(mu**(1.0_dp-power)-1.0_dp)/(1.0_dp-power)
if(abs(power-2.0_dp)<1.0e-7_dp)then
kappa=log(mu)+(2.0_dp-power)*log(mu)**2/2.0_dp
else
kappa=(mu**(2.0_dp-power)-1.0_dp)/(2.0_dp-power)
end if
m1=exp((y*theta-kappa)/phi)
m2=1.0_dp/mu
dev=tweedie_dev(y,mu,power)
m3=exp(-dev/(2.0_dp*phi))/y
im=0_c_int
if(present(method))im=int(method,c_int)
if(im<1_c_int.or.im>3_c_int)then
 im=1_c_int
 if(m2<m1)im=2_c_int
 if(m3<merge(m1,m2,im==1_c_int))im=3_c_int
end if
select case(im)
case(1)
ya(1)=y
phia(1)=phi
scale=m1
case(2)
ya(1)=y/mu
phia(1)=phi/mu**(2.0_dp-power)
scale=m2
case default
ya(1)=1.0_dp
phia(1)=phi/y**(2.0_dp-power)
scale=m3
end select
mua(1)=1.0_dp
pp=power
n=1_c_int
vb=merge(1_c_int,0_c_int,verb)
pdf=1_c_int
call twcomputation_loop(n,pp,phia,ya,mua,vb,pdf,fv,es,re,its)
d=real(fv(1),dp)*scale
call set_details(int(es(1)),real(re(1),dp),int(its(1)))
contains
subroutine set_details(e,r,it)
integer,intent(in)::e,it
real(dp),intent(in)::r
if(present(exitstatus))exitstatus=e
if(present(relerr))relerr=r
if(present(regions))regions=it
end subroutine
end function

function ptweedie_inversion(q,mu,phi,power,verbose,exitstatus,relerr,regions,igexact) result(pval)
real(dp),intent(in)::q,mu,phi,power
logical,intent(in),optional::verbose,igexact
integer,intent(out),optional::exitstatus,regions
real(dp),intent(out),optional::relerr
real(dp)::pval
real(c_double)::qa(1),mua(1),phia(1),fv(1),re(1),pp
integer(c_int)::es(1),its(1),n,vb,pdf
logical::verb,ig
if(phi<=0.0_dp.or.mu<=0.0_dp.or.power<1.0_dp)then
pval=ieee_value(0.0_dp,ieee_quiet_nan)
call set_details(-1,0.0_dp,0)
return
end if
verb=.false.
if(present(verbose))verb=verbose
ig=.true.
if(present(igexact))ig=igexact
if(power==1.0_dp)then
pval=ppois(q/phi,lambda=mu/phi)
call set_details(0,0.0_dp,0)
return
end if
if(power==2.0_dp)then
pval=pgamma(q,shape=1.0_dp/phi,rate=1.0_dp/(phi*mu))
call set_details(0,0.0_dp,0)
return
end if
if(power==3.0_dp.and.ig)then
pval=pinvgauss_tweedie(q,mu,phi)
call set_details(0,0.0_dp,0)
return
end if
if(q<0.0_dp)then
pval=0.0_dp
call set_details(0,0.0_dp,0)
return
end if
if(q==0.0_dp)then
 if(power<2.0_dp)then
 pval=exp(-mu**(2.0_dp-power)/(phi*(2.0_dp-power)))
 else
 pval=0.0_dp
 end if
 call set_details(0,0.0_dp,0)
 return
end if
qa(1)=q
mua(1)=mu
phia(1)=phi
pp=power
n=1_c_int
vb=merge(1_c_int,0_c_int,verb)
pdf=0_c_int
call twcomputation_loop(n,pp,phia,qa,mua,vb,pdf,fv,es,re,its)
pval=min(1.0_dp,max(0.0_dp,real(fv(1),dp)))
call set_details(int(es(1)),real(re(1),dp),int(its(1)))
contains
subroutine set_details(e,r,it)
integer,intent(in)::e,it
real(dp),intent(in)::r
if(present(exitstatus))exitstatus=e
if(present(relerr))relerr=r
if(present(regions))regions=it
end subroutine
end function

end module tweedie_inversion_mod
