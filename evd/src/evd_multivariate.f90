! SPDX-License-Identifier: GPL-3.0-only
module evd_multivariate
use r_compat, only: dp, runif1, rexp
use evd_transform, only: gev_to_exp_measure, exp_measure_to_gev
implicit none
private
public :: pmvlog, dmvlog, rmvlog, amvlog
public :: pmvalog, dmvalog, rmvalog, amvalog
public :: log_positive_stable
contains

pure function exponent_alog(x, dep, asy) result(v)
real(dp), intent(in) :: x(:), dep(:), asy(:,:)
real(dp) :: v,s
integer :: b,j
v=0.0_dp
do b=1,size(dep)
 s=0.0_dp
 do j=1,size(x)
  if(asy(b,j)>0.0_dp .and. x(j)>0.0_dp) s=s+(asy(b,j)*x(j))**(1.0_dp/dep(b))
 end do
 if(s>0.0_dp) v=v+s**dep(b)
end do
end function

pure function pmvlog(q,dep,mar,lower_tail) result(p)
real(dp),intent(in)::q(:),dep,mar(:,:)
logical,intent(in),optional::lower_tail
real(dp)::p,x(size(q)),v
logical::lt
integer::j
lt=.true.
if(present(lower_tail))lt=lower_tail
do j=1,size(q)
x(j)=gev_to_exp_measure(q(j),mar(min(j,size(mar,1)),1),mar(min(j,size(mar,1)),2),mar(min(j,size(mar,1)),3))
end do
v=sum(x**(1.0_dp/dep))**dep
p=exp(-v)
if(.not.lt) p=upper_tail_inclusion_log(q,dep,mar)
end function

pure recursive function upper_tail_inclusion_log(q,dep,mar) result(p)
real(dp),intent(in)::q(:),dep,mar(:,:)
real(dp)::p,tmp(size(q))
integer::mask,j,d,bits
d=size(q)
p=0.0_dp
do mask=0,2**d-1
 tmp=q
 bits=0
 do j=1,d
  if(btest(mask,j-1))then
  tmp(j)=huge(1.0_dp)
  bits=bits+1
  end if
 end do
 p=p+(-1.0_dp)**(d-bits)*pmvlog(tmp,dep,mar,.true.)
end do
end function

pure function pmvalog(q,dep,asy,mar,lower_tail) result(p)
real(dp),intent(in)::q(:),dep(:),asy(:,:),mar(:,:)
logical,intent(in),optional::lower_tail
real(dp)::p,x(size(q)),v
integer::j
logical::lt
lt=.true.
if(present(lower_tail))lt=lower_tail
do j=1,size(q)
x(j)=gev_to_exp_measure(q(j),mar(min(j,size(mar,1)),1),mar(min(j,size(mar,1)),2),mar(min(j,size(mar,1)),3))
end do
v=exponent_alog(x,dep,asy)
p=exp(-v)
if(.not.lt) p=upper_tail_inclusion_alog(q,dep,asy,mar)
end function

pure recursive function upper_tail_inclusion_alog(q,dep,asy,mar) result(p)
real(dp),intent(in)::q(:),dep(:),asy(:,:),mar(:,:)
real(dp)::p,tmp(size(q))
integer::mask,j,d,bits
d=size(q)
p=0.0_dp
do mask=0,2**d-1
 tmp=q
 bits=0
 do j=1,d
  if(btest(mask,j-1))then
  tmp(j)=huge(1.0_dp)
  bits=bits+1
  end if
 end do
 p=p+(-1.0_dp)**(d-bits)*pmvalog(tmp,dep,asy,mar,.true.)
end do
end function

pure function block_derivative(mask,x,dep,asy) result(vd)
integer,intent(in)::mask
real(dp),intent(in)::x(:),dep(:),asy(:,:)
real(dp)::vd,s,coef,prod,r
integer::b,j,k,m
logical::ok
m=popcnt(mask)
vd=0.0_dp
if(m==0)return
do b=1,size(dep)
 ok=.true.
 s=0.0_dp
 do j=1,size(x)
  if(asy(b,j)>0.0_dp .and. x(j)>0.0_dp) s=s+(asy(b,j)*x(j))**(1.0_dp/dep(b))
  if(btest(mask,j-1).and.(asy(b,j)<=0.0_dp.or.x(j)<=0.0_dp))ok=.false.
 end do
 if(.not.ok.or.s<=0.0_dp)cycle
 coef=1.0_dp
 do k=0,m-1
 coef=coef*(dep(b)-real(k,dp))
 end do
 if(abs(coef)<=tiny(1.0_dp))cycle
 r=1.0_dp/dep(b)
 coef=coef*r**m*s**(dep(b)-real(m,dp))
 prod=1.0_dp
 do j=1,size(x)
  if(btest(mask,j-1))prod=prod*asy(b,j)**r*x(j)**(r-1.0_dp)
 end do
 vd=vd+coef*prod
end do
end function

recursive subroutine partition_sum(pos,d,rgs,maxb,x,dep,asy,total)
integer,intent(in)::pos,d,maxb
integer,intent(inout)::rgs(:)
real(dp),intent(in)::x(:),dep(:),asy(:,:)
real(dp),intent(inout)::total
integer::b,nb,j,mask
real(dp)::prod,vd
if(pos>d)then
 nb=maxval(rgs(1:d))
 prod=1.0_dp
 do b=1,nb
  mask=0
  do j=1,d
   if(rgs(j)==b)mask=ibset(mask,j-1)
  end do
  vd=block_derivative(mask,x,dep,asy)
  prod=prod*(-vd)
 end do
 total=total+prod
 return
end if
do b=1,maxb+1
 rgs(pos)=b
 call partition_sum(pos+1,d,rgs,max(maxb,b),x,dep,asy,total)
end do
end subroutine

function dmvalog(xobs,dep,asy,mar,log_) result(dns)
real(dp),intent(in)::xobs(:),dep(:),asy(:,:),mar(:,:)
logical,intent(in),optional::log_
real(dp)::dns,x(size(xobs)),v,psum,jac
integer::j,d,rgs(size(xobs))
logical::lg
lg=.false.
if(present(log_))lg=log_
d=size(xobs)
do j=1,d
 x(j)=gev_to_exp_measure(xobs(j),mar(min(j,size(mar,1)),1),mar(min(j,size(mar,1)),2),mar(min(j,size(mar,1)),3))
 if(x(j)<=0.0_dp)then
 dns=merge(-huge(1.0_dp),0.0_dp,lg)
 return
 end if
end do
v=exponent_alog(x,dep,asy)
rgs=1
psum=0.0_dp
rgs(1)=1
call partition_sum(2,d,rgs,1,x,dep,asy,psum)
if(mod(d,2)==1) psum=-psum
jac=0.0_dp
do j=1,d
 jac=jac+(1.0_dp+mar(min(j,size(mar,1)),3))*log(x(j))-log(mar(min(j,size(mar,1)),2))
end do
if(psum<=0.0_dp)then
dns=merge(-huge(1.0_dp),0.0_dp,lg)
else
dns=log(psum)-v+jac
if(.not.lg)dns=exp(dns)
end if
end function

function dmvlog(xobs,dep,mar,log_) result(dns)
real(dp),intent(in)::xobs(:),dep,mar(:,:)
logical,intent(in),optional::log_
real(dp)::dns
integer::d,j
real(dp),allocatable::deps(:),asy(:,:)
d=size(xobs)
allocate(deps(1),asy(1,d))
deps=dep
asy=1.0_dp
dns=dmvalog(xobs,deps,asy,mar,log_)
end function

pure function amvlog(x,dep) result(a)
real(dp),intent(in)::x(:),dep
real(dp)::a,s
s=sum(x)
if(s<=0.0_dp)then
a=0.0_dp
else
a=sum((x/s)**(1.0_dp/dep))**dep
end if
end function
pure function amvalog(x,dep,asy) result(a)
real(dp),intent(in)::x(:),dep(:),asy(:,:)
real(dp)::a,s,xx(size(x))
s=sum(x)
if(s<=0.0_dp)then
a=0.0_dp
return
end if
xx=x/s
a=exponent_alog(xx,dep,asy)
end function

function log_positive_stable(alpha) result(s)
! Returns the logarithm of the positive stable variate used by Stephenson's algorithms.
real(dp),intent(in)::alpha
real(dp)::s,t,u,w,a,pi
if(alpha==1.0_dp)then
s=0.0_dp
return
end if
pi=acos(-1.0_dp)
t=1.0_dp-alpha
u=pi*runif1()
w=log(-log(max(runif1(),tiny(1.0_dp))))
a=log(sin(t*u))+(alpha/t)*log(sin(alpha*u))-(1.0_dp/t)*log(sin(u))
s=(t/alpha)*(a-w)
end function

function rmvlog(n,d,dep,mar) result(sim)
integer,intent(in)::n,d
real(dp),intent(in)::dep,mar(:,:)
real(dp)::sim(n,d),s,e
integer::i,j
do i=1,n
 s=log_positive_stable(dep)
 do j=1,d
  e=-log(max(runif1(),tiny(1.0_dp)))
  sim(i,j)=exp(dep*(s-log(e)))
  sim(i,j)=exp_measure_to_gev(1.0_dp/sim(i,j),mar(min(j,size(mar,1)),1),mar(min(j,size(mar,1)),2),mar(min(j,size(mar,1)),3))
 end do
end do
end function

function rmvalog(n,d,dep,asy,mar) result(sim)
integer,intent(in)::n,d
real(dp),intent(in)::dep(:),asy(:,:),mar(:,:)
real(dp)::sim(n,d),gev(size(dep),d),s,e,z
integer::i,j,b
sim=0.0_dp
do i=1,n
 gev=0.0_dp
 do b=1,size(dep)
  if(dep(b)/=1.0_dp)then
  s=log_positive_stable(dep(b))
  else
  s=0.0_dp
  end if
  do j=1,d
   if(asy(b,j)>0.0_dp)then
    e=-log(max(runif1(),tiny(1.0_dp)))
    gev(b,j)=asy(b,j)*exp(dep(b)*(s-log(e)))
   end if
  end do
 end do
 do j=1,d
  z=maxval(gev(:,j))
  sim(i,j)=exp_measure_to_gev(1.0_dp/z,mar(min(j,size(mar,1)),1),mar(min(j,size(mar,1)),2),mar(min(j,size(mar,1)),3))
 end do
end do
end function

end module evd_multivariate
