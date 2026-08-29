! SPDX-License-Identifier: GPL-3.0-only
module evd_bivariate
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_negative_inf
use r_compat, only: dp, normal_cdf, dnorm, pbeta, dbeta
use evd_transform, only: gev_to_exp_measure
implicit none
private
public :: pbvlog,pbvalog,pbvhr,pbvneglog,pbvaneglog,pbvbilog,pbvnegbilog,pbvct,pbvamix
public :: dbvlog,dbvalog,dbvhr,dbvneglog,dbvaneglog,dbvbilog,dbvnegbilog,dbvct,dbvamix
public :: abvlog,abvalog,abvhr,abvneglog,abvaneglog,abvbilog,abvnegbilog,abvct,abvamix
public :: hbvlog,hbvalog,hbvhr,hbvneglog,hbvaneglog,hbvbilog,hbvnegbilog,hbvct,hbvamix
public :: ccbvlog,ccbvalog,ccbvhr,ccbvneglog,ccbvaneglog,ccbvbilog,ccbvnegbilog,ccbvct,ccbvamix
public :: bilog_gamma, negbilog_gamma
contains

pure function finish_prob(v,x1,x2,lower_tail) result(p)
real(dp),intent(in)::v,x1,x2
logical,intent(in),optional::lower_tail
real(dp)::p
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
p=exp(-v)
if(.not.lt)p=1.0_dp-exp(-x1)-exp(-x2)+p
p=max(0.0_dp,min(1.0_dp,p))
end function

pure subroutine tx(q1,q2,mar1,mar2,x1,x2)
real(dp),intent(in)::q1,q2,mar1(3),mar2(3)
real(dp),intent(out)::x1,x2
x1=gev_to_exp_measure(q1,mar1(1),mar1(2),mar1(3))
x2=gev_to_exp_measure(q2,mar2(1),mar2(2),mar2(3))
end subroutine

pure function pbvlog(q1,q2,dep,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,dep,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,v
call tx(q1,q2,mar1,mar2,x1,x2)
v=(x1**(1.0_dp/dep)+x2**(1.0_dp/dep))**dep
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvalog(q1,q2,dep,asy,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,dep,asy(2),mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,v
call tx(q1,q2,mar1,mar2,x1,x2)
v=((asy(1)*x1)**(1.0_dp/dep)+(asy(2)*x2)**(1.0_dp/dep))**dep + &
  (1.0_dp-asy(1))*x1+(1.0_dp-asy(2))*x2
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvhr(q1,q2,dep,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,dep,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,v
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
 v=x1+x2
else
 v=x1*normal_cdf(1.0_dp/dep+dep*log(x1/x2)/2.0_dp)+ &
   x2*normal_cdf(1.0_dp/dep+dep*log(x2/x1)/2.0_dp)
end if
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvneglog(q1,q2,dep,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,dep,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,v,z
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
z=0.0_dp
else
z=(x1**(-dep)+x2**(-dep))**(-1.0_dp/dep)
end if
v=x1+x2-z
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvaneglog(q1,q2,dep,asy,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,dep,asy(2),mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,v,z
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp.or.any(asy<=0.0_dp))then
z=0.0_dp
else
z=((asy(1)*x1)**(-dep)+(asy(2)*x2)**(-dep))**(-1.0_dp/dep)
end if
v=x1+x2-z
p=finish_prob(v,x1,x2,lower_tail)
end function

pure function bilog_gamma(x1,x2,alpha,beta) result(g)
real(dp),intent(in)::x1,x2,alpha,beta
real(dp)::g,lo,hi,mid,flo,fm
integer::i
if(x1==0.0_dp)then
g=0.0_dp
return
end if
if(x2==0.0_dp)then
g=1.0_dp
return
end if
lo=0.0_dp
hi=1.0_dp
flo=(1.0_dp-alpha)*x1
 do i=1,100
  mid=0.5_dp*(lo+hi)
  fm=(1.0_dp-alpha)*x1*(1.0_dp-mid)**beta-(1.0_dp-beta)*x2*mid**alpha
  if(fm==0.0_dp.or.hi-lo<sqrt(epsilon(1.0_dp)))exit
  if(sign(1.0_dp,flo)/=sign(1.0_dp,fm))then
  hi=mid
  else
  lo=mid
  flo=fm
  end if
 end do
g=0.5_dp*(lo+hi)
end function
pure function negbilog_gamma(x1,x2,alpha,beta) result(g)
real(dp),intent(in)::x1,x2,alpha,beta
real(dp)::g,lo,hi,mid,flo,fm
integer::i
if(x1==0.0_dp)then
g=1.0_dp
return
end if
if(x2==0.0_dp)then
g=0.0_dp
return
end if
lo=0.0_dp
hi=1.0_dp
flo=-(1.0_dp+beta)*x2
 do i=1,100
  mid=0.5_dp*(lo+hi)
  fm=(1.0_dp+alpha)*x1*mid**alpha-(1.0_dp+beta)*x2*(1.0_dp-mid)**beta
  if(fm==0.0_dp.or.hi-lo<sqrt(epsilon(1.0_dp)))exit
  if(sign(1.0_dp,flo)/=sign(1.0_dp,fm))then
  hi=mid
  else
  lo=mid
  flo=fm
  end if
 end do
g=0.5_dp*(lo+hi)
end function

pure function pbvbilog(q1,q2,alpha,beta,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,g,v
call tx(q1,q2,mar1,mar2,x1,x2)
g=bilog_gamma(x1,x2,alpha,beta)
v=x1*g**(1.0_dp-alpha)+x2*(1.0_dp-g)**(1.0_dp-beta)
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvnegbilog(q1,q2,alpha,beta,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,g,v
call tx(q1,q2,mar1,mar2,x1,x2)
g=negbilog_gamma(x1,x2,alpha,beta)
v=x1+x2-x1*g**(1.0_dp+alpha)-x2*(1.0_dp-g)**(1.0_dp+beta)
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvct(q1,q2,alpha,beta,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,u,v
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1+x2==0.0_dp)then
v=0.0_dp
else
 u=alpha*x2/(alpha*x2+beta*x1)
 v=x2*pbeta(u,alpha,beta+1.0_dp)+x1*(1.0_dp-pbeta(u,alpha+1.0_dp,beta))
end if
p=finish_prob(v,x1,x2,lower_tail)
end function
pure function pbvamix(q1,q2,alpha,beta,mar1,mar2,lower_tail) result(p)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::lower_tail
real(dp)::p,x1,x2,s,v
call tx(q1,q2,mar1,mar2,x1,x2)
s=x1+x2
if(s==0.0_dp)then
v=0.0_dp
else
v=s-(alpha+beta)*x1+alpha*x1*x1/s+beta*x1**3/s**2
end if
p=finish_prob(v,x1,x2,lower_tail)
end function

pure function jacobian_log(x1,x2,mar1,mar2) result(j)
real(dp),intent(in)::x1,x2,mar1(3),mar2(3)
real(dp)::j
j=(1.0_dp+mar1(3))*log(x1)+(1.0_dp+mar2(3))*log(x2)-log(mar1(2)*mar2(2))
end function

pure function dbvlog(q1,q2,dep,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,dep,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,idep,z,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 idep=1.0_dp/dep
 z=(x1**idep+x2**idep)**dep
 ld=(idep+mar1(3))*log(x1)+(idep+mar2(3))*log(x2)-log(mar1(2)*mar2(2)) + &
    (1.0_dp-2.0_dp*idep)*log(z)+log(idep-1.0_dp+z)-z
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvalog(q1,q2,dep,asy,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,dep,asy(2),mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,idep,z,v,j,e1,e2,e3,e4,e5,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp.or.any(asy<=0.0_dp).or.any(asy>=1.0_dp))then
 ! Boundary asymmetric cases are distributions with singular components; the ordinary 2-D density is not represented by this formula.
 ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 idep=1.0_dp/dep
 z=((asy(1)*x1)**idep+(asy(2)*x2)**idep)**dep
 v=z+(1.0_dp-asy(1))*x1+(1.0_dp-asy(2))*x2
 j=jacobian_log(x1,x2,mar1,mar2)
 e1=log(1.0_dp-asy(1))+log(1.0_dp-asy(2))
 e2=log(1.0_dp-asy(1))+idep*log(asy(2))+(idep-1.0_dp)*log(x2)
 e3=log(1.0_dp-asy(2))+idep*log(asy(1))+(idep-1.0_dp)*log(x1)
 e4=(1.0_dp-idep)*log(z)+log(exp(e2)+exp(e3))
 e5=idep*log(asy(1))+idep*log(asy(2))+(idep-1.0_dp)*(log(x1)+log(x2)) + &
    (1.0_dp-2.0_dp*idep)*log(z)+log(idep-1.0_dp+z)
 ld=log(exp(e1)+exp(e4)+exp(e5))-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvhr(q1,q2,dep,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,dep,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,a,b,v,e,j,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 a=1.0_dp/dep+dep*log(x1/x2)/2.0_dp
 b=1.0_dp/dep+dep*log(x2/x1)/2.0_dp
 v=x1*normal_cdf(a)+x2*normal_cdf(b)
 e=x1*normal_cdf(a)*x2*normal_cdf(b)+dep*x1*dnorm(a)/2.0_dp
 j=mar1(3)*log(x1)+mar2(3)*log(x2)-log(mar1(2)*mar2(2))
 ld=log(e)+j-v
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvneglog(q1,q2,dep,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,dep,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,z,v,e1,e2,j,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 z=(x1**(-dep)+x2**(-dep))**(-1.0_dp/dep)
 v=x1+x2-z
 j=jacobian_log(x1,x2,mar1,mar2)
 e1=(1.0_dp+dep)*log(z)+log(exp((-dep-1.0_dp)*log(x1))+exp((-dep-1.0_dp)*log(x2)))
 e2=(-dep-1.0_dp)*(log(x1)+log(x2))+(1.0_dp+2.0_dp*dep)*log(z)+log(1.0_dp+dep+z)
 ld=log(1.0_dp-exp(e1)+exp(e2))-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvaneglog(q1,q2,dep,asy,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,dep,asy(2),mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,z,v,e1,e2,e3,e4,j,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp.or.any(asy<=0.0_dp))then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 z=((asy(1)*x1)**(-dep)+(asy(2)*x2)**(-dep))**(-1.0_dp/dep)
 v=x1+x2-z
 j=jacobian_log(x1,x2,mar1,mar2)
 e1=-dep*log(asy(1))+(-dep-1.0_dp)*log(x1)
 e2=-dep*log(asy(2))+(-dep-1.0_dp)*log(x2)
 e3=(1.0_dp+dep)*log(z)+log(exp(e1)+exp(e2))
 e4=-dep*(log(asy(1))+log(asy(2)))+(-dep-1.0_dp)*(log(x1)+log(x2)) + &
    (1.0_dp+2.0_dp*dep)*log(z)+log(1.0_dp+dep+z)
 ld=log(1.0_dp-exp(e3)+exp(e4))-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvbilog(q1,q2,alpha,beta,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,g,v,j,e1,e2,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 g=bilog_gamma(x1,x2,alpha,beta)
 v=x1*g**(1.0_dp-alpha)+x2*(1.0_dp-g)**(1.0_dp-beta)
 j=jacobian_log(x1,x2,mar1,mar2)
 e1=g**(1.0_dp-alpha)*(1.0_dp-g)**(1.0_dp-beta)
 e2=(1.0_dp-alpha)*beta*(1.0_dp-g)**(beta-1.0_dp)*x1 + &
    (1.0_dp-beta)*alpha*g**(alpha-1.0_dp)*x2
 ld=log(e1+(1.0_dp-alpha)*(1.0_dp-beta)/e2)-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvnegbilog(q1,q2,alpha,beta,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,g,v,j,e1,e2,e3,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 g=negbilog_gamma(x1,x2,alpha,beta)
 v=x1+x2-x1*g**(1.0_dp+alpha)-x2*(1.0_dp-g)**(1.0_dp+beta)
 j=jacobian_log(x1,x2,mar1,mar2)
 e1=(1.0_dp-g**(1.0_dp+alpha))*(1.0_dp-(1.0_dp-g)**(1.0_dp+beta))
 e2=(1.0_dp+alpha)*(1.0_dp+beta)*g**alpha*(1.0_dp-g)**beta
 e3=(1.0_dp+alpha)*alpha*g**(alpha-1.0_dp)*x1 + &
    (1.0_dp+beta)*beta*(1.0_dp-g)**(beta-1.0_dp)*x2
 ld=log(e1+e2/e3)-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvct(q1,q2,alpha,beta,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,u,v,j,c1,e1,e2,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 u=alpha*x2/(alpha*x2+beta*x1)
 v=x2*pbeta(u,alpha,beta+1.0_dp)+x1*(1.0_dp-pbeta(u,alpha+1.0_dp,beta))
 j=jacobian_log(x1,x2,mar1,mar2)
 c1=alpha*beta/(alpha+beta+1.0_dp)
 e1=pbeta(u,alpha,beta+1.0_dp)*(1.0_dp-pbeta(u,alpha+1.0_dp,beta))
 e2=dbeta(u,alpha+1.0_dp,beta+1.0_dp)/(alpha*x2+beta*x1)
 ld=log(e1+c1*e2)-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure function dbvamix(q1,q2,alpha,beta,mar1,mar2,log_) result(d)
real(dp),intent(in)::q1,q2,alpha,beta,mar1(3),mar2(3)
logical,intent(in),optional::log_
real(dp)::d,x1,x2,s,v,j,u1,u2,v1,v2,v12,ld
logical::lg
lg=.false.
if(present(log_))lg=log_
call tx(q1,q2,mar1,mar2,x1,x2)
if(x1<=0.0_dp.or.x2<=0.0_dp)then
ld=ieee_value(0.0_dp,ieee_negative_inf)
else
 s=x1+x2
 v=s-(alpha+beta)*x1+alpha*x1*x1/s+beta*x1**3/s**2
 j=jacobian_log(x1,x2,mar1,mar2)
 u1=x1/s
 u2=x2/s
 v1=1.0_dp-alpha*u2**2-beta*(3.0_dp*u2**2-2.0_dp*u2**3)
 v2=1.0_dp-alpha*u1**2-2.0_dp*beta*u1**3
 v12=(-2.0_dp*alpha*u1*u2-6.0_dp*beta*u1**2*u2)/s
 ld=log(v1*v2-v12)-v+j
end if
if(lg)then
d=ld
else
d=exp(ld)
end if
end function

pure elemental function abvlog(x,dep) result(a)
real(dp),intent(in)::x,dep
real(dp)::a
a=(x**(1.0_dp/dep)+(1.0_dp-x)**(1.0_dp/dep))**dep
end function
pure elemental function abvalog(x,dep,asy1,asy2) result(a)
real(dp),intent(in)::x,dep,asy1,asy2
real(dp)::a
a=((asy1*x)**(1.0_dp/dep)+(asy2*(1.0_dp-x))**(1.0_dp/dep))**dep+(1.0_dp-asy1)*x+(1.0_dp-asy2)*(1.0_dp-x)
end function
pure elemental function abvhr(x,dep) result(a)
real(dp),intent(in)::x,dep
real(dp)::a
if(x<=0.0_dp.or.x>=1.0_dp)then
a=1.0_dp
else
a=x*normal_cdf(1.0_dp/dep+dep*log(x/(1.0_dp-x))/2.0_dp) + &
 (1.0_dp-x)*normal_cdf(1.0_dp/dep+dep*log((1.0_dp-x)/x)/2.0_dp)
 end if
end function
pure elemental function abvneglog(x,dep) result(a)
real(dp),intent(in)::x,dep
real(dp)::a
if(x<=0.0_dp.or.x>=1.0_dp)then
a=1.0_dp
else
a=1.0_dp-(x**(-dep)+(1.0_dp-x)**(-dep))**(-1.0_dp/dep)
end if
end function
pure elemental function abvaneglog(x,dep,asy1,asy2) result(a)
real(dp),intent(in)::x,dep,asy1,asy2
real(dp)::a
if(x<=0.0_dp.or.x>=1.0_dp.or.asy1==0.0_dp.or.asy2==0.0_dp)then
a=1.0_dp
else
a=1.0_dp-((asy1*x)**(-dep)+(asy2*(1.0_dp-x))**(-dep))**(-1.0_dp/dep)
end if
end function
pure function abvbilog(x,alpha,beta) result(a)
real(dp),intent(in)::x,alpha,beta
real(dp)::a,g
g=bilog_gamma(x,1.0_dp-x,alpha,beta)
a=x*g**(1.0_dp-alpha)+(1.0_dp-x)*(1.0_dp-g)**(1.0_dp-beta)
end function
pure function abvnegbilog(x,alpha,beta) result(a)
real(dp),intent(in)::x,alpha,beta
real(dp)::a,g
g=negbilog_gamma(x,1.0_dp-x,alpha,beta)
a=1.0_dp-x*g**(1.0_dp+alpha)-(1.0_dp-x)*(1.0_dp-g)**(1.0_dp+beta)
end function
pure function abvct(x,alpha,beta) result(a)
real(dp),intent(in)::x,alpha,beta
real(dp)::a,u
if(x<=0.0_dp.or.x>=1.0_dp)then
a=1.0_dp
return
end if
u=alpha*(1.0_dp-x)/(alpha*(1.0_dp-x)+beta*x)
a=(1.0_dp-x)*pbeta(u,alpha,beta+1.0_dp)+x*(1.0_dp-pbeta(u,alpha+1.0_dp,beta))
end function
pure elemental function abvamix(x,alpha,beta) result(a)
real(dp),intent(in)::x,alpha,beta
real(dp)::a
a=1.0_dp-(alpha+beta)*x+alpha*x*x+beta*x**3
end function

pure elemental function hbvlog(x,dep,half) result(h)
real(dp),intent(in)::x,dep
logical,intent(in),optional::half
real(dp)::h,idep
if(x<=0.0_dp.or.x>=1.0_dp)then
h=0.0_dp
return
end if
idep=1.0_dp/dep
h=(idep-1.0_dp)*(x*(1.0_dp-x))**(-1.0_dp-idep) * &
 (x**(-idep)+(1.0_dp-x)**(-idep))**(dep-2.0_dp)
 if(present(half))then
 if(half)h=h/2
 end if
end function
pure elemental function hbvalog(x,dep,asy1,asy2,half) result(h)
real(dp),intent(in)::x,dep,asy1,asy2
logical,intent(in),optional::half
real(dp)::h,idep
if(x<=0.0_dp.or.x>=1.0_dp.or.asy1==0.0_dp.or.asy2==0.0_dp)then
h=0.0_dp
return
end if
idep=1.0_dp/dep
h=(idep-1.0_dp)*(asy1*asy2)**idep*(x*(1.0_dp-x))**(-1.0_dp-idep)* &
 ((asy1/x)**idep+(asy2/(1.0_dp-x))**idep)**(dep-2.0_dp)
 if(present(half))then
 if(half)h=h/2
 end if
end function
pure elemental function hbvhr(x,dep,half) result(h)
real(dp),intent(in)::x,dep
logical,intent(in),optional::half
real(dp)::h,z
if(x<=0.0_dp.or.x>=1.0_dp)then
h=0.0_dp
return
end if
z=1.0_dp/dep+dep*log(x/(1.0_dp-x))/2.0_dp
h=dep*dnorm(z)/(2.0_dp*x*(1.0_dp-x)**2)
if(present(half))then
if(half)h=h/2
end if
end function
pure elemental function hbvneglog(x,dep,half) result(h)
real(dp),intent(in)::x,dep
logical,intent(in),optional::half
real(dp)::h
if(x<=0.0_dp.or.x>=1.0_dp)then
h=0.0_dp
return
end if
h=(1.0_dp+dep)*(x*(1.0_dp-x))**(dep-1.0_dp)*(x**dep+(1.0_dp-x)**dep)**(-1.0_dp/dep-2.0_dp)
if(present(half))then
if(half)h=h/2
end if
end function
pure elemental function hbvaneglog(x,dep,asy1,asy2,half) result(h)
real(dp),intent(in)::x,dep,asy1,asy2
logical,intent(in),optional::half
real(dp)::h
if(x<=0.0_dp.or.x>=1.0_dp.or.asy1==0.0_dp.or.asy2==0.0_dp)then
h=0.0_dp
return
end if
h=(1.0_dp+dep)*(asy1*asy2)**(-dep)*(x*(1.0_dp-x))**(dep-1.0_dp)* &
 ((x/asy1)**dep+((1.0_dp-x)/asy2)**dep)**(-1.0_dp/dep-2.0_dp)
if(present(half))then
if(half)h=h/2
end if
end function
pure function hbvbilog(x,alpha,beta,half) result(h)
real(dp),intent(in)::x,alpha,beta
logical,intent(in),optional::half
real(dp)::h,g,den
if(x<=0.0_dp.or.x>=1.0_dp)then
h=0.0_dp
return
end if
g=bilog_gamma(1.0_dp-x,x,alpha,beta)
den=(1.0_dp-alpha)*beta*(1.0_dp-g)**(beta-1.0_dp)*(1.0_dp-x)+(1.0_dp-beta)*alpha*g**(alpha-1.0_dp)*x
h=(1.0_dp-alpha)*(1.0_dp-beta)/(x*(1.0_dp-x)*den)
if(present(half))then
if(half)h=h/2
end if
end function
pure function hbvnegbilog(x,alpha,beta,half) result(h)
real(dp),intent(in)::x,alpha,beta
logical,intent(in),optional::half
real(dp)::h,g,den
if(x<=0.0_dp.or.x>=1.0_dp)then
h=0.0_dp
return
end if
g=negbilog_gamma(1.0_dp-x,x,alpha,beta)
den=(1.0_dp+alpha)*alpha*g**(alpha-1.0_dp)*(1.0_dp-x)+(1.0_dp+beta)*beta*(1.0_dp-g)**(beta-1.0_dp)*x
h=(1.0_dp+alpha)*(1.0_dp+beta)*g**alpha*(1.0_dp-g)**beta/(x*(1.0_dp-x)*den)
if(present(half))then
if(half)h=h/2
end if
end function
pure function hbvct(x,alpha,beta,half) result(h)
real(dp),intent(in)::x,alpha,beta
logical,intent(in),optional::half
real(dp)::h,u,c1
if(x<=0.0_dp.or.x>=1.0_dp)then
h=0.0_dp
return
end if
u=alpha*x/(alpha*x+beta*(1.0_dp-x))
c1=alpha*beta/(alpha+beta+1.0_dp)
h=dbeta(u,alpha+1.0_dp,beta+1.0_dp)/(alpha*x*x*(1.0_dp-x)+beta*x*(1.0_dp-x)**2)*c1
if(present(half))then
if(half)h=h/2
end if
end function
pure elemental function hbvamix(x,alpha,beta,half) result(h)
real(dp),intent(in)::x,alpha,beta
logical,intent(in),optional::half
real(dp)::h
h=2.0_dp*alpha+6.0_dp*beta*(1.0_dp-x)
if(present(half))then
if(half)h=h/2
end if
end function

pure function ccbvlog(m1,m2,oldm1,dep) result(f)
real(dp),intent(in)::m1,m2,oldm1,dep
real(dp)::f,t1,t2,u,v,idep
t1=-log(m1)
t2=-log(m2)
idep=1.0_dp/dep
u=t1**idep+t2**idep
v=u**dep
f=exp(-v)/m2*t2**(idep-1.0_dp)*u**(dep-1.0_dp)-oldm1
end function
pure function ccbvalog(m1,m2,oldm1,dep,asy1,asy2) result(f)
real(dp),intent(in)::m1,m2,oldm1,dep,asy1,asy2
real(dp)::f,t1,t2,u,v,idep
t1=-log(m1)
t2=-log(m2)
idep=1.0_dp/dep
u=(asy1*t1)**idep+(asy2*t2)**idep
v=(1.0_dp-asy1)*t1+(1.0_dp-asy2)*t2+u**dep
f=exp(-v)/m2*(1.0_dp-asy2+asy2**idep*t2**(idep-1.0_dp)*u**(dep-1.0_dp))-oldm1
end function
pure function ccbvhr(m1,m2,oldm1,dep) result(f)
real(dp),intent(in)::m1,m2,oldm1,dep
real(dp)::f,t1,t2,v,z
t1=-log(m1)
t2=-log(m2)
z=1.0_dp/dep+(log(t2)-log(t1))*dep/2.0_dp
v=t2*normal_cdf(z)+t1*normal_cdf(1.0_dp/dep+(log(t1)-log(t2))*dep/2.0_dp)
f=normal_cdf(z)*exp(-v)/m2-oldm1
end function
pure function ccbvneglog(m1,m2,oldm1,dep) result(f)
real(dp),intent(in)::m1,m2,oldm1,dep
real(dp)::f,t1,t2,v,idep
t1=-log(m1)
t2=-log(m2)
idep=1.0_dp/dep
v=(t2**(-dep)+t1**(-dep))**(-idep)
f=exp(v)*m1*(1.0_dp-(1.0_dp+(t2/t1)**dep)**(-1.0_dp-idep))-oldm1
end function
pure function ccbvaneglog(m1,m2,oldm1,dep,asy1,asy2) result(f)
real(dp),intent(in)::m1,m2,oldm1,dep,asy1,asy2
real(dp)::f,t1,t2,v,idep
t1=-log(m1)
t2=-log(m2)
idep=1.0_dp/dep
v=(asy1*t2)**(-dep)+(asy2*t1)**(-dep)
f=exp(v**(-idep))*m1*(1.0_dp-asy1**(-dep)*t2**(-dep-1.0_dp)*v**(-idep-1.0_dp))-oldm1
end function
pure function ccbvbilog(m1,m2,oldm1,alpha,beta) result(f)
real(dp),intent(in)::m1,m2,oldm1,alpha,beta
real(dp)::f,t1,t2,g,v
t1=-log(m1)
t2=-log(m2)
g=bilog_gamma(t1,t2,alpha,beta)
v=t1*g**(1.0_dp-alpha)+t2*(1.0_dp-g)**(1.0_dp-beta)
f=exp(-v)/m2*(1.0_dp-g)**(1.0_dp-beta)-oldm1
end function
pure function ccbvnegbilog(m1,m2,oldm1,alpha,beta) result(f)
real(dp),intent(in)::m1,m2,oldm1,alpha,beta
real(dp)::f,t1,t2,g,v
t1=-log(m1)
t2=-log(m2)
g=negbilog_gamma(t1,t2,alpha,beta)
v=-t1-t2+t1*g**(1.0_dp+alpha)+t2*(1.0_dp-g)**(1.0_dp+beta)
f=exp(v)/m2*(1.0_dp-(1.0_dp-g)**(1.0_dp+beta))-oldm1
end function
pure function ccbvct(m1,m2,oldm1,alpha,beta) result(f)
real(dp),intent(in)::m1,m2,oldm1,alpha,beta
real(dp)::f,t1,t2,u,v
t1=-log(m1)
t2=-log(m2)
u=alpha*t2/(alpha*t2+beta*t1)
v=t1*(1.0_dp-pbeta(u,alpha+1.0_dp,beta))+t2*pbeta(u,alpha,beta+1.0_dp)
f=exp(-v)/m2*pbeta(u,alpha,beta+1.0_dp)-oldm1
end function
pure function ccbvamix(m1,m2,oldm1,alpha,beta) result(f)
real(dp),intent(in)::m1,m2,oldm1,alpha,beta
real(dp)::f,t1,t2,u,v,v2
t1=-log(m1)
t2=-log(m2)
u=t1/(t1+t2)
v=t1+t2-t1*((alpha+beta)-alpha*u-beta*u*u)
v2=1.0_dp-alpha*u*u-2.0_dp*beta*u**3
f=exp(-v)/m2*v2-oldm1
end function

end module evd_bivariate
