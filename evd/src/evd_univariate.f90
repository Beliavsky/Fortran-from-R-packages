! SPDX-License-Identifier: GPL-3.0-only
module evd_univariate
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
use r_compat, only: dp, runif1, rexp, rbeta, r_lgamma
implicit none
private
public :: dgev, pgev, qgev, rgev
public :: dgumbel, pgumbel, qgumbel, rgumbel
public :: dfrechet, pfrechet, qfrechet, rfrechet
public :: drweibull, prweibull, qrweibull, rrweibull
public :: dnweibull, pnweibull, qnweibull, rnweibull
public :: dgpd, pgpd, qgpd, rgpd
public :: dgumbelx, pgumbelx, qgumbelx, rgumbelx
public :: pextreme_from_cdf, dextreme_from_pdf_cdf, qextreme_to_base_prob
public :: porder_from_cdf, dorder_from_pdf_cdf, rorder_base_prob, rextreme_base_prob

contains

pure elemental function dgev(x, loc, scale, shape, log_) result(d)
real(dp), intent(in) :: x, loc, scale, shape
logical, intent(in), optional :: log_
real(dp) :: d, z, t, ld
logical :: lg
lg = .false.
if (present(log_)) lg = log_
if (scale <= 0.0_dp) then
   d = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
z = (x-loc)/scale
if (shape == 0.0_dp) then
   ld = -log(scale) - z - exp(-z)
else
   t = 1.0_dp + shape*z
   if (t <= 0.0_dp) then
      ld = ieee_value(0.0_dp, ieee_negative_inf)
   else
      ld = -log(scale) - t**(-1.0_dp/shape) - (1.0_dp/shape + 1.0_dp)*log(t)
   end if
end if
if (lg) then
d=ld
else
d=exp(ld)
end if
end function dgev

pure elemental function pgev(q, loc, scale, shape, lower_tail) result(p)
real(dp), intent(in) :: q, loc, scale, shape
logical, intent(in), optional :: lower_tail
real(dp) :: p, z, t
logical :: lt
lt=.true.
if(present(lower_tail)) lt=lower_tail
if(scale <= 0.0_dp) then
   p=ieee_value(0.0_dp,ieee_quiet_nan)
   return
end if
z=(q-loc)/scale
if(shape==0.0_dp) then
   p=exp(-exp(-z))
else
   t=max(1.0_dp+shape*z,0.0_dp)
   if(t==0.0_dp) then
      if(shape>0.0_dp) then
      p=0.0_dp
      else
      p=1.0_dp
      end if
   else
      p=exp(-t**(-1.0_dp/shape))
   end if
end if
if(.not.lt) p=1.0_dp-p
end function pgev

pure elemental function qgev(p, loc, scale, shape, lower_tail) result(q)
real(dp), intent(in) :: p, loc, scale, shape
logical, intent(in), optional :: lower_tail
real(dp) :: q, pp
logical :: lt
lt=.true.
if(present(lower_tail)) lt=lower_tail
pp=p
if(.not.lt) pp=1.0_dp-pp
if(scale < 0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
   q=ieee_value(0.0_dp,ieee_quiet_nan)
   return
end if
if(shape==0.0_dp) then
   q=loc-scale*log(-log(pp))
else
   q=loc+scale*((-log(pp))**(-shape)-1.0_dp)/shape
end if
end function qgev

function rgev(n, loc, scale, shape) result(x)
integer, intent(in) :: n
real(dp), intent(in) :: loc, scale, shape
real(dp) :: x(n), e(n)
e = rexp(n, 1.0_dp)
if(shape==0.0_dp) then
   x=loc-scale*log(e)
else
   x=loc+scale*(e**(-shape)-1.0_dp)/shape
end if
end function rgev

pure elemental function dgumbel(x, loc, scale, log_) result(d)
real(dp),intent(in)::x,loc,scale
logical,intent(in),optional::log_
real(dp)::d
d=dgev(x,loc,scale,0.0_dp,log_)
end function
pure elemental function pgumbel(q,loc,scale,lower_tail) result(p)
real(dp),intent(in)::q,loc,scale
logical,intent(in),optional::lower_tail
real(dp)::p
p=pgev(q,loc,scale,0.0_dp,lower_tail)
end function
pure elemental function qgumbel(p,loc,scale,lower_tail) result(q)
real(dp),intent(in)::p,loc,scale
logical,intent(in),optional::lower_tail
real(dp)::q
q=qgev(p,loc,scale,0.0_dp,lower_tail)
end function
function rgumbel(n,loc,scale) result(x)
integer,intent(in)::n
real(dp),intent(in)::loc,scale
real(dp)::x(n)
x=rgev(n,loc,scale,0.0_dp)
end function

pure elemental function dfrechet(x,loc,scale,shape,log_) result(d)
real(dp),intent(in)::x,loc,scale,shape
logical,intent(in),optional::log_
real(dp)::d,z,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
if(scale<=0.0_dp.or.shape<=0.0_dp)then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z=(x-loc)/scale
if(z<=0.0_dp)then
 ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 ld=log(shape/scale)-(1.0_dp+shape)*log(z)-z**(-shape)
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function
pure elemental function pfrechet(q,loc,scale,shape,lower_tail) result(p)
real(dp),intent(in)::q,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::p,z
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
if(scale<=0.0_dp.or.shape<=0.0_dp)then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z=max((q-loc)/scale,0.0_dp)
if(z==0.0_dp)then
p=0.0_dp
else
p=exp(-z**(-shape))
end if
if(.not.lt)p=1.0_dp-p
end function
pure elemental function qfrechet(p,loc,scale,shape,lower_tail) result(q)
real(dp),intent(in)::p,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::q,pp
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
pp=p
if(.not.lt)pp=1.0_dp-pp
if(scale<0.0_dp.or.shape<=0.0_dp.or.pp<=0.0_dp.or.pp>=1.0_dp)then
q=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
q=loc+scale*(-log(pp))**(-1.0_dp/shape)
end function
function rfrechet(n,loc,scale,shape) result(x)
integer,intent(in)::n
real(dp),intent(in)::loc,scale,shape
real(dp)::x(n),e(n)
e=rexp(n,1.0_dp)
x=loc+scale*e**(-1.0_dp/shape)
end function

pure elemental function drweibull(x,loc,scale,shape,log_) result(d)
real(dp),intent(in)::x,loc,scale,shape
logical,intent(in),optional::log_
real(dp)::d,z,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
if(scale<=0.0_dp.or.shape<=0.0_dp)then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z=(x-loc)/scale
if(z>=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
ld=log(shape/scale)+(shape-1.0_dp)*log(-z)-(-z)**shape
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function
pure elemental function prweibull(q,loc,scale,shape,lower_tail) result(p)
real(dp),intent(in)::q,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::p,z
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
if(scale<=0.0_dp.or.shape<=0.0_dp)then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z=min((q-loc)/scale,0.0_dp)
p=exp(-(-z)**shape)
if(.not.lt)p=1.0_dp-p
end function
pure elemental function qrweibull(p,loc,scale,shape,lower_tail) result(q)
real(dp),intent(in)::p,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::q,pp
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
pp=p
if(.not.lt)pp=1.0_dp-pp
if(scale<0.0_dp.or.shape<=0.0_dp.or.pp<=0.0_dp.or.pp>=1.0_dp)then
q=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
q=loc-scale*(-log(pp))**(1.0_dp/shape)
end function
function rrweibull(n,loc,scale,shape) result(x)
integer,intent(in)::n
real(dp),intent(in)::loc,scale,shape
real(dp)::x(n),e(n)
e=rexp(n,1.0_dp)
x=loc-scale*e**(1.0_dp/shape)
end function

pure elemental function dnweibull(x,loc,scale,shape,log_) result(d)
real(dp),intent(in)::x,loc,scale,shape
logical,intent(in),optional::log_
real(dp)::d
d=drweibull(x,loc,scale,shape,log_)
end function
pure elemental function pnweibull(q,loc,scale,shape,lower_tail) result(p)
real(dp),intent(in)::q,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::p
p=prweibull(q,loc,scale,shape,lower_tail)
end function
pure elemental function qnweibull(p,loc,scale,shape,lower_tail) result(q)
real(dp),intent(in)::p,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::q
q=qrweibull(p,loc,scale,shape,lower_tail)
end function
function rnweibull(n,loc,scale,shape) result(x)
integer,intent(in)::n
real(dp),intent(in)::loc,scale,shape
real(dp)::x(n)
x=rrweibull(n,loc,scale,shape)
end function

pure elemental function dgpd(x,loc,scale,shape,log_) result(d)
real(dp),intent(in)::x,loc,scale,shape
logical,intent(in),optional::log_
real(dp)::d,z,t,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
if(scale<=0.0_dp)then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z=(x-loc)/scale
t=1.0_dp+shape*z
if(z<=0.0_dp.or.t<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else if(shape==0.0_dp)then
ld=-log(scale)-z
else
ld=-log(scale)-(1.0_dp/shape+1.0_dp)*log(t)
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function
pure elemental function pgpd(q,loc,scale,shape,lower_tail) result(p)
real(dp),intent(in)::q,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::p,z,t
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
if(scale<=0.0_dp)then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z=max(q-loc,0.0_dp)/scale
if(shape==0.0_dp)then
p=1.0_dp-exp(-z)
else
t=max(1.0_dp+shape*z,0.0_dp)
if(t==0.0_dp)then
p=1.0_dp
else
p=1.0_dp-t**(-1.0_dp/shape)
end if
end if
if(.not.lt)p=1.0_dp-p
end function
pure elemental function qgpd(p,loc,scale,shape,lower_tail) result(q)
real(dp),intent(in)::p,loc,scale,shape
logical,intent(in),optional::lower_tail
real(dp)::q,pp
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
pp=p
if(lt)pp=1.0_dp-pp
if(scale<0.0_dp.or.pp<=0.0_dp.or.pp>=1.0_dp)then
q=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(shape==0.0_dp)then
q=loc-scale*log(pp)
else
q=loc+scale*(pp**(-shape)-1.0_dp)/shape
end if
end function
function rgpd(n,loc,scale,shape) result(x)
integer,intent(in)::n
real(dp),intent(in)::loc,scale,shape
real(dp)::x(n),u
integer::i
do i=1,n
 u=runif1()
 if(shape==0.0_dp)then
 x(i)=loc-scale*log(u)
 else
 x(i)=loc+scale*(u**(-shape)-1.0_dp)/shape
 end if
end do
end function

pure elemental function pgumbelx(q,loc1,scale1,loc2,scale2,lower_tail) result(p)
real(dp),intent(in)::q,loc1,scale1,loc2,scale2
logical,intent(in),optional::lower_tail
real(dp)::p,z1,z2
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
if(scale1<=0.0_dp.or.scale2<=0.0_dp.or.loc1>loc2)then
p=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z1=(q-loc1)/scale1
z2=(q-loc2)/scale2
p=exp(-exp(-z1)-exp(-z2))
if(.not.lt)p=1.0_dp-p
end function
pure elemental function dgumbelx(x,loc1,scale1,loc2,scale2,log_) result(d)
real(dp),intent(in)::x,loc1,scale1,loc2,scale2
logical,intent(in),optional::log_
real(dp)::d,z1,z2
logical::lg
lg=.false.
if(present(log_))lg=log_
if(scale1<=0.0_dp.or.scale2<=0.0_dp.or.loc1>loc2)then
d=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
z1=(x-loc1)/scale1
z2=(x-loc2)/scale2
d=exp(-exp(-z1)-z2-exp(-z2))/scale2+exp(-exp(-z2)-z1-exp(-z1))/scale1
if(lg)d=log(d)
end function
function rgumbelx(n,loc1,scale1,loc2,scale2) result(x)
integer,intent(in)::n
real(dp),intent(in)::loc1,scale1,loc2,scale2
real(dp)::x(n),a(n),b(n)
a=rgumbel(n,loc1,scale1)
b=rgumbel(n,loc2,scale2)
x=max(a,b)
end function
function qgumbelx(p,lower,upper,loc1,scale1,loc2,scale2,lower_tail,tol) result(q)
real(dp),intent(in)::p,lower,upper,loc1,scale1,loc2,scale2
logical,intent(in),optional::lower_tail
real(dp),intent(in),optional::tol
real(dp)::q,lo,hi,mid,target,eps
integer::i
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
target=p
if(.not.lt)target=1.0_dp-p
eps=sqrt(epsilon(1.0_dp))
if(present(tol))eps=tol
lo=lower
hi=upper
do i=1,200
 mid=0.5_dp*(lo+hi)
 if(pgumbelx(mid,loc1,scale1,loc2,scale2)<target)then
 lo=mid
 else
 hi=mid
 end if
 if(abs(hi-lo)<=eps*max(1.0_dp,abs(mid)))exit
end do
q=0.5_dp*(lo+hi)
end function

pure elemental function pextreme_from_cdf(f,mlen,largest,lower_tail) result(p)
real(dp),intent(in)::f
integer,intent(in)::mlen
logical,intent(in),optional::largest,lower_tail
real(dp)::p,g
logical::la,lt
la=.true.
lt=.true.
if(present(largest))la=largest
if(present(lower_tail))lt=lower_tail
g=f
if(.not.la)g=1.0_dp-g
p=g**mlen
if(la.neqv.lt)p=1.0_dp-p
end function

pure elemental function dextreme_from_pdf_cdf(pdf,cdf,mlen,largest,log_) result(d)
real(dp),intent(in)::pdf,cdf
integer,intent(in)::mlen
logical,intent(in),optional::largest,log_
real(dp)::d,g,ld
logical::la,lg
la=.true.
lg=.false.
if(present(largest))la=largest
if(present(log_))lg=log_
g=cdf
if(.not.la)g=1.0_dp-g
if(pdf<=0.0_dp.or.g<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
ld=log(real(mlen,dp))+log(pdf)+real(mlen-1,dp)*log(g)
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure elemental function qextreme_to_base_prob(p,mlen,largest,lower_tail) result(pb)
real(dp),intent(in)::p
integer,intent(in)::mlen
logical,intent(in),optional::largest,lower_tail
real(dp)::pb,pp
logical::la,lt
la=.true.
lt=.true.
if(present(largest))la=largest
if(present(lower_tail))lt=lower_tail
pp=p
if(.not.lt)pp=1.0_dp-pp
if(la)then
pb=pp**(1.0_dp/real(mlen,dp))
else
pb=1.0_dp-(1.0_dp-pp)**(1.0_dp/real(mlen,dp))
end if
end function

pure function porder_from_cdf(f,mlen,j,largest,lower_tail) result(p)
real(dp),intent(in)::f
integer,intent(in)::mlen,j
logical,intent(in),optional::largest,lower_tail
real(dp)::p,term,lf
integer::k,k0,k1
logical::la,lt
la=.true.
lt=.true.
if(present(largest))la=largest
if(present(lower_tail))lt=lower_tail
p=0.0_dp
if(la)then
k0=mlen+1-j
k1=mlen
else
k0=0
k1=j-1
end if
do k=k0,k1
 lf=r_lgamma(real(mlen+1,dp))-r_lgamma(real(k+1,dp))-r_lgamma(real(mlen-k+1,dp))
 if(f==0.0_dp.and.k>0)cycle
 if(f==1.0_dp.and.mlen-k>0)cycle
 term=exp(lf)
 if(k>0)term=term*f**k
 if(mlen-k>0)term=term*(1.0_dp-f)**(mlen-k)
 p=p+term
end do
if(la.neqv.lt)p=1.0_dp-p
end function

pure function dorder_from_pdf_cdf(pdf,cdf,mlen,j,largest,log_) result(d)
real(dp),intent(in)::pdf,cdf
integer,intent(in)::mlen,j
logical,intent(in),optional::largest,log_
real(dp)::d,ld
integer::jj
logical::la,lg
la=.true.
lg=.false.
if(present(largest))la=largest
if(present(log_))lg=log_
jj=j
if(.not.la)jj=mlen+1-jj
if(pdf<=0.0_dp.or.cdf<=0.0_dp.or.cdf>=1.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 ld=r_lgamma(real(mlen+1,dp))-r_lgamma(real(jj,dp))-r_lgamma(real(mlen-jj+1,dp)) + &
    log(pdf)+real(mlen-jj,dp)*log(cdf)+real(jj-1,dp)*log(1.0_dp-cdf)
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

function rextreme_base_prob(n,mlen,largest) result(p)
integer,intent(in)::n,mlen
logical,intent(in),optional::largest
real(dp)::p(n)
logical::la
la=.true.
if(present(largest))la=largest
if(la)then
p=rbeta(n,real(mlen,dp),1.0_dp)
else
p=rbeta(n,1.0_dp,real(mlen,dp))
end if
end function
function rorder_base_prob(n,mlen,j,largest) result(p)
integer,intent(in)::n,mlen,j
logical,intent(in),optional::largest
real(dp)::p(n)
integer::jj
logical::la
la=.true.
if(present(largest))la=largest
jj=j
if(.not.la)jj=mlen+1-j
p=rbeta(n,real(mlen+1-jj,dp),real(jj,dp))
end function

end module evd_univariate
