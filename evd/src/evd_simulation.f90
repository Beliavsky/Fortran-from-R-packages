! SPDX-License-Identifier: GPL-3.0-only
module evd_simulation
use r_compat, only: dp, runif1
use evd_transform, only: exp_measure_to_gev
use evd_bivariate, only: ccbvlog,ccbvalog,ccbvhr,ccbvneglog,ccbvaneglog,ccbvbilog,ccbvene=>ccbvnegbilog,ccbvc=>ccbvct,ccbva=>ccbvamix
implicit none
private
public :: rbvlog,rbvalog,rbvhr,rbvneglog,rbvaneglog,rbvbilog,rbvnegbilog,rbvct,rbvamix,evmc
contains

function rbvlog(n,dep,mar1,mar2) result(sim)
! Stephenson (2003), Algorithm 1.1, as used by upstream evd.
integer,intent(in)::n
real(dp),intent(in)::dep,mar1(3),mar2(3)
real(dp)::sim(n,2),u,z,e1,e2
integer::i
do i=1,n
 u=runif1()
 if(runif1()<dep)then
 e1=-log(runif1())
 e2=-log(runif1())
 z=e1+e2
 else
 z=-log(runif1())
 end if
 sim(i,1)=1.0_dp/(z*u**dep)
 sim(i,2)=1.0_dp/(z*(1.0_dp-u)**dep)
 sim(i,1)=exp_measure_to_gev(1.0_dp/sim(i,1),mar1(1),mar1(2),mar1(3))
 sim(i,2)=exp_measure_to_gev(1.0_dp/sim(i,2),mar2(1),mar2(2),mar2(3))
end do
end function

function rbvalog(n,dep,asy,mar1,mar2) result(sim)
! Stephenson (2003), Algorithm 1.2, as used by upstream evd.
integer,intent(in)::n
real(dp),intent(in)::dep,asy(2),mar1(3),mar2(3)
real(dp)::sim(n,2),u,z,v11,v22,v12,v21,e1,e2
integer::i
if(dep==1.0_dp.or.any(asy==0.0_dp))then
 do i=1,n
  sim(i,1)=exp_measure_to_gev(-log(runif1()),mar1(1),mar1(2),mar1(3))
  sim(i,2)=exp_measure_to_gev(-log(runif1()),mar2(1),mar2(2),mar2(3))
 end do
 return
end if
do i=1,n
 v11=(1.0_dp-asy(1))/(-log(runif1()))
 v22=(1.0_dp-asy(2))/(-log(runif1()))
 u=runif1()
 if(runif1()<dep)then
 e1=-log(runif1())
 e2=-log(runif1())
 z=e1+e2
 else
 z=-log(runif1())
 end if
 v12=asy(1)/(z*u**dep)
 v21=asy(2)/(z*(1.0_dp-u)**dep)
 sim(i,1)=max(v11,v12)
 sim(i,2)=max(v22,v21)
 sim(i,1)=exp_measure_to_gev(1.0_dp/sim(i,1),mar1(1),mar1(2),mar1(3))
 sim(i,2)=exp_measure_to_gev(1.0_dp/sim(i,2),mar2(1),mar2(2),mar2(3))
end do
end function

pure function cond_value(model,m1,m2,old,dep,asy,alpha,beta) result(f)
character(len=*),intent(in)::model
real(dp),intent(in)::m1,m2,old,dep,asy(2),alpha,beta
real(dp)::f
select case(trim(model))
case('log');f=ccbvlog(m1,m2,old,dep)
case('alog');f=ccbvalog(m1,m2,old,dep,asy(1),asy(2))
case('hr');f=ccbvhr(m1,m2,old,dep)
case('neglog');f=ccbvneglog(m1,m2,old,dep)
case('aneglog');f=ccbvaneglog(m1,m2,old,dep,asy(1),asy(2))
case('bilog');f=ccbvbilog(m1,m2,old,alpha,beta)
case('negbilog');f=ccbvene(m1,m2,old,alpha,beta)
case('ct');f=ccbvc(m1,m2,old,alpha,beta)
case('amix');f=ccbva(m1,m2,old,alpha,beta)
case default;f=0.0_dp
end select
end function

function conditional_uniform_pair(old1,m2,model,dep,asy,alpha,beta) result(m1)
real(dp),intent(in)::old1,m2,dep,asy(2),alpha,beta
character(len=*),intent(in)::model
real(dp)::m1,lo,hi,mid,flo,fm,eps
integer::j
eps=sqrt(epsilon(1.0_dp))
lo=eps
hi=1.0_dp-eps
flo=cond_value(model,lo,m2,old1,dep,asy,alpha,beta)
do j=1,60
 mid=0.5_dp*(lo+hi)
 fm=cond_value(model,mid,m2,old1,dep,asy,alpha,beta)
 if(abs(fm)<eps.or.hi-lo<eps)exit
 if(sign(1.0_dp,flo)/=sign(1.0_dp,fm))then
 hi=mid
 else
 lo=mid
 flo=fm
 end if
end do
m1=0.5_dp*(lo+hi)
end function

function rbv_conditional(n,model,dep,asy,alpha,beta,mar1,mar2) result(sim)
integer,intent(in)::n
character(len=*),intent(in)::model
real(dp),intent(in)::dep,asy(2),alpha,beta,mar1(3),mar2(3)
real(dp)::sim(n,2),u1,u2
integer::i
do i=1,n
 u1=runif1()
 u2=runif1()
 u1=conditional_uniform_pair(u1,u2,model,dep,asy,alpha,beta)
 sim(i,1)=exp_measure_to_gev(-log(u1),mar1(1),mar1(2),mar1(3))
 sim(i,2)=exp_measure_to_gev(-log(u2),mar2(1),mar2(2),mar2(3))
end do
end function

function rbvhr(n,dep,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::dep,mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'hr',dep,[1.0_dp,1.0_dp],1.0_dp,1.0_dp,mar1,mar2)
end function
function rbvneglog(n,dep,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::dep,mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'neglog',dep,[1.0_dp,1.0_dp],1.0_dp,1.0_dp,mar1,mar2)
end function
function rbvaneglog(n,dep,asy,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::dep,asy(2),mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'aneglog',dep,asy,1.0_dp,1.0_dp,mar1,mar2)
end function
function rbvbilog(n,alpha,beta,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::alpha,beta,mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'bilog',1.0_dp,[1.0_dp,1.0_dp],alpha,beta,mar1,mar2)
end function
function rbvnegbilog(n,alpha,beta,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::alpha,beta,mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'negbilog',1.0_dp,[1.0_dp,1.0_dp],alpha,beta,mar1,mar2)
end function
function rbvct(n,alpha,beta,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::alpha,beta,mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'ct',1.0_dp,[1.0_dp,1.0_dp],alpha,beta,mar1,mar2)
end function
function rbvamix(n,alpha,beta,mar1,mar2) result(sim)
integer,intent(in)::n
real(dp),intent(in)::alpha,beta,mar1(3),mar2(3)
real(dp)::sim(n,2)
sim=rbv_conditional(n,'amix',1.0_dp,[1.0_dp,1.0_dp],alpha,beta,mar1,mar2)
end function

function evmc(n,model,dep,asy,alpha,beta,margin) result(x)
integer,intent(in)::n
character(len=*),intent(in)::model,margin
real(dp),intent(in)::dep,asy(2),alpha,beta
real(dp)::x(n),u1,u2
integer::i
x(1)=runif1()
do i=2,n
 u1=runif1()
 u2=x(i-1)
 x(i)=conditional_uniform_pair(u1,u2,model,dep,asy,alpha,beta)
end do
select case(trim(margin))
case('uniform')
case('frechet');x=-1.0_dp/log(x)
case('rweibull');x=log(x)
case('gumbel');x=-log(-log(x))
end select
end function
end module evd_simulation
