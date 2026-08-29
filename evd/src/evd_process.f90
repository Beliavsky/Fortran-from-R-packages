! SPDX-License-Identifier: GPL-3.0-only
module evd_process
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use r_compat, only: dp
use evd_univariate, only: rfrechet
implicit none
private
public :: cluster_result_t, clusters, exi, marma, mar, mma

type :: cluster_result_t
 integer :: nclusters=0
 integer, allocatable :: label(:)
 logical, allocatable :: is_start(:),is_end(:)
 real(dp) :: acs=0.0_dp
 real(dp), allocatable :: maxima(:)
end type
contains

function clusters(data,u,r,ulow,rlow) result(out)
real(dp),intent(in)::data(:),u(:),ulow(:)
integer,intent(in)::r,rlow
type(cluster_result_t)::out
integer::n,i,j,rr,clind,shigh,shigh2,k
logical::incl
logical,allocatable::high(:),high2(:)
real(dp),allocatable::mx(:)
n=size(data)
allocate(out%label(n),out%is_start(n),out%is_end(n),high(n),high2(n))
out%label=0
out%is_start=.false.
out%is_end=.false.
high=data>u
high2=data>ulow
incl=.false.
clind=0
 do i=1,n
  if(high(i).and.incl)out%label(i)=clind
  if(high(i).and..not.incl)then
   incl=.true.
   clind=clind+1
   out%is_start(i)=.true.
   out%label(i)=clind
  end if
  if(.not.high(i).and.incl)then
   rr=min(r,n-i+1)
   shigh=0
   do j=i,i+rr-1
   if(high(j))shigh=shigh+1
   end do
   rr=min(rlow,n-i+1)
   shigh2=0
   do j=i,i+rr-1
   if(high2(j))shigh2=shigh2+1
   end do
   if(shigh==0.or.shigh2==0)then
    incl=.false.
    if(i>1)out%is_end(i-1)=.true.
   else
    out%label(i)=clind
   end if
  end if
 end do
if(incl)out%is_end(n)=.true.
out%nclusters=clind
if(clind>0)then
 out%acs=real(count(high),dp)/real(clind,dp)
 allocate(mx(clind))
 mx=-huge(1.0_dp)
 do i=1,n
  k=out%label(i)
  if(k>0)mx(k)=max(mx(k),data(i))
 end do
 out%maxima=mx
else
 out%acs=0.0_dp
 allocate(out%maxima(0))
end if
end function

function exi(data,u,r,ulow,rlow) result(theta)
real(dp),intent(in)::data(:),u(:),ulow(:)
integer,intent(in)::r,rlow
real(dp)::theta,den
integer,allocatable::idx(:),gap(:)
integer::i,nx,k
type(cluster_result_t)::cl
if(r>=1)then
 cl=clusters(data,u,r,ulow,rlow)
 if(cl%acs>0.0_dp)then
 theta=1.0_dp/cl%acs
 else
 theta=ieee_value(0.0_dp,ieee_quiet_nan)
 end if
 return
end if
nx=count(data>u)
if(nx==0)then
theta=ieee_value(0.0_dp,ieee_quiet_nan)
return
end if
if(nx==1)then
theta=1.0_dp
return
end if
allocate(idx(nx))
k=0
do i=1,size(data)
if(data(i)>u(i))then
k=k+1
idx(k)=i
end if
end do
allocate(gap(nx-1))
gap=idx(2:)-idx(:nx-1)
if(maxval(gap)>2)then
 den=log(real(nx-1,dp))+log(sum(real((gap-1)*(gap-2),dp)))
 theta=min(1.0_dp,exp(log(2.0_dp)+2.0_dp*log(sum(real(gap-1,dp)))-den))
else
 den=log(real(nx-1,dp))+log(sum(real(gap*gap,dp)))
 theta=min(1.0_dp,exp(log(2.0_dp)+2.0_dp*log(sum(real(gap,dp)))-den))
end if
end function

function marma(n,psi,theta,init,n_start,innov) result(x)
integer,intent(in)::n,n_start
real(dp),intent(in)::psi(:),theta(:),init(:),innov(:)
real(dp)::x(n),work(size(init)+n+n_start-size(init)),th0(size(theta)+1),v
integer::p,q,i,j,ix
p=size(psi)
q=size(theta)
work=0.0_dp
if(p>0)work(1:p)=init
th0(1)=1.0_dp
if(q>0)th0(2:)=theta
do i=1,n+n_start-p
 ix=i+p
 v=0.0_dp
 do j=1,p
 v=max(v,psi(j)*work(ix-j))
 end do
 do j=0,q
 v=max(v,th0(j+1)*innov(i+q-j))
 end do
 work(ix)=v
end do
x=work(n_start+1:n_start+n)
end function

function mar(n,psi,init,n_start) result(x)
integer,intent(in)::n,n_start
real(dp),intent(in)::psi(:),init(:)
real(dp)::x(n),innov(n+n_start)
innov=rfrechet(n+n_start,0.0_dp,1.0_dp,1.0_dp)
x=marma(n,psi,[real(dp)::],init,n_start,innov)
end function

function mma(n,theta) result(x)
integer,intent(in)::n
real(dp),intent(in)::theta(:)
real(dp)::x(n),innov(n+size(theta))
innov=rfrechet(n+size(theta),0.0_dp,1.0_dp,1.0_dp)
x=marma(n,[real(dp)::],theta,[real(dp)::],0,innov)
end function
end module evd_process
