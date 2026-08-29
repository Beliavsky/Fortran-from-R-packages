! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_ate
use r_compat, only: dp, normal_cdf, dnorm
implicit none
private
public :: ate_moments, ate_gradient, ate_marginal_effects
public :: ATE_BAL,ATE_BAL_SAMPLE,ATE_ATT,ATE_LINEAR,ATE_LOGIT,ATE_PROBIT
integer,parameter :: ATE_BAL=1,ATE_BAL_SAMPLE=2,ATE_ATT=3
integer,parameter :: ATE_LINEAR=1,ATE_LOGIT=2,ATE_PROBIT=3
contains

pure subroutine ate_moments(theta,x,k,mom_type,family,gt,w,pop_mom)
! Computational kernel of .momentFctATE. x layout follows upstream:
! y | treatment design (k columns, including baseline column) | balance covariates.
real(dp),intent(in)::theta(:),x(:,:)
integer,intent(in)::k,mom_type,family
real(dp),allocatable,intent(out)::gt(:,:)
real(dp),intent(in),optional::w(:,:),pop_mom(:)
real(dp),allocatable::z(:,:),covar(:,:),gt1(:,:),gt2(:,:),gt3(:,:),extra(:,:)
real(dp),allocatable::tetz(:),tetb(:),eta(:),mu(:),e(:),target(:)
integer::n,l,nh,i,j,col,qbase
n=size(x,1)
nh=size(x,2)-(k+1)
if(present(w))then
 allocate(z(n,k+size(w,2)))
 z(:,1:k)=x(:,2:k+1)
 z(:,k+1:)=w
else
 allocate(z(n,k))
 z=x(:,2:k+1)
end if
l=size(z,2)
tetz=theta(1:l)
tetb=theta(size(theta)-k+2:size(theta))
eta=matmul(z,tetz)
allocate(mu(n))
select case(family)
case(ATE_LOGIT);mu=1.0_dp/(1.0_dp+exp(-eta))
case(ATE_PROBIT)
do i=1,n
mu(i)=normal_cdf(eta(i))
end do
case default;mu=eta
end select
e=x(:,1)-mu
gt1=z*spread(e,2,l)
gt2=x(:,3:k+1)-spread(tetb,1,n)
if(nh>0)then
 covar=x(:,k+2:)
 allocate(gt3(n,nh*(k-1)))
 col=0
 do i=1,k-1
 do j=1,nh
 col=col+1
 gt3(:,col)=gt2(:,i)*covar(:,j)
 end do
 end do
else
allocate(covar(n,0),gt3(n,0))
end if
qbase=l+(k-1)+nh*(k-1)
allocate(gt(n,qbase))
gt(:,1:l)=gt1
gt(:,l+1:l+k-1)=gt2
if(size(gt3,2)>0)gt(:,l+k:)=gt3
if(nh<=0)return
if(present(pop_mom))then
 allocate(extra(n,nh))
 extra=covar-spread(pop_mom,1,n)
 call append_cols(gt,extra)
else if(mom_type==ATE_BAL_SAMPLE)then
 allocate(target(nh))
 target=sum(covar,dim=1)/real(n,dp)
 allocate(extra(n,nh))
 extra=covar-spread(target,1,n)
 call append_cols(gt,extra)
else if(mom_type==ATE_ATT)then
 ! Upstream ATT target is the treated-group mean for each treatment column.
 allocate(extra(n,nh))
 extra=0.0_dp
 do j=1,nh
   block
    real(dp)::den,targ
    den=sum(x(:,3))
    if(den>0)then
    targ=sum(x(:,3)*covar(:,j))/den
    else
    targ=0
    end if
    extra(:,j)=covar(:,j)-targ
   end block
 end do
 call append_cols(gt,extra)
end if
end subroutine ate_moments

pure subroutine ate_gradient(theta,x,k,mom_type,family,g,w,pop_mom,pt)
! Weighted Jacobian of the ATE moment vector; mirrors .DmomentFctATE.
real(dp),intent(in)::theta(:),x(:,:)
integer,intent(in)::k,mom_type,family
real(dp),allocatable,intent(out)::g(:,:)
real(dp),intent(in),optional::w(:,:),pop_mom(:),pt(:)
real(dp),allocatable::z(:,:),covar(:,:),p(:),eta(:),tau(:)
integer::n,l,nh,q,nt,i,j,row
n=size(x,1)
nh=size(x,2)-(k+1)
if(present(w))then
allocate(z(n,k+size(w,2)))
z(:,1:k)=x(:,2:k+1)
z(:,k+1:)=w
else
allocate(z(n,k))
z=x(:,2:k+1)
end if
l=size(z,2)
nt=size(theta)
eta=matmul(z,theta(1:l))
allocate(tau(n))
select case(family)
case(ATE_LOGIT);tau=exp(-eta)/(1.0_dp+exp(-eta))**2
case(ATE_PROBIT)
do i=1,n
tau(i)=dnorm(eta(i))
end do
case default;tau=1.0_dp
end select
allocate(p(n))
if(present(pt))then
p=pt
else
p=1.0_dp/real(n,dp)
end if
q=l+k-1+nh*(k-1)
if(mom_type/=ATE_BAL.or.present(pop_mom))q=q+nh
allocate(g(q,nt))
g=0.0_dp
do i=1,l
do j=1,l
g(i,j)=-sum(p*z(:,i)*tau*z(:,j))
end do
end do
do i=1,k-1
g(l+i,l+i)=-sum(p)
end do
if(nh>0)then
 covar=x(:,k+2:)
 row=l+k
 do i=1,k-1
 do j=1,nh
   g(row,l+i)=-sum(p*covar(:,j))
   row=row+1
 end do
 end do
end if
end subroutine ate_gradient

subroutine ate_marginal_effects(coef,vcov,k,family,estimate,se)
real(dp),intent(in)::coef(:),vcov(:,:)
integer,intent(in)::k,family
real(dp),allocatable,intent(out)::estimate(:),se(:)
real(dp)::p0,d0,p1,d1,a(2),v2(2,2)
integer::i
allocate(estimate(k),se(k))
call link_values(coef(1),family,p0,d0)
estimate(1)=p0
se(1)=abs(d0)*sqrt(max(vcov(1,1),0.0_dp))
do i=1,k-1
 call link_values(coef(1)+coef(i+1),family,p1,d1)
 estimate(i+1)=p1-p0
 a=[d1-d0,d1]
 v2=vcov([1,i+1],[1,i+1])
 se(i+1)=sqrt(max(dot_product(a,matmul(v2,a)),0.0_dp))
end do
end subroutine ate_marginal_effects

pure subroutine link_values(eta,family,mu,dmu)
real(dp),intent(in)::eta
integer,intent(in)::family
real(dp),intent(out)::mu,dmu
select case(family)
case(ATE_LOGIT)
mu=1/(1+exp(-eta))
dmu=mu*(1-mu)
case(ATE_PROBIT)
mu=normal_cdf(eta)
dmu=dnorm(eta)
case default
mu=eta
dmu=1
end select
end subroutine link_values

pure subroutine append_cols(a,b)
real(dp),allocatable,intent(inout)::a(:,:)
real(dp),intent(in)::b(:,:)
real(dp),allocatable::tmp(:,:)
allocate(tmp(size(a,1),size(a,2)+size(b,2)))
tmp(:,1:size(a,2))=a
tmp(:,size(a,2)+1:)=b
call move_alloc(tmp,a)
end subroutine append_cols
end module gmm_ate
