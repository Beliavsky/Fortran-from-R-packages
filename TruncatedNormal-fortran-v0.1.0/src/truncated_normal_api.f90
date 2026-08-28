! R-style computational wrappers for TruncatedNormal.
! SPDX-License-Identifier: GPL-3.0-only
module truncated_normal_api
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
 use r_mod, only: dp, solve_real
 use truncated_normal_math, only: lnNpr,phinv,norminvp,trandn,qtnorm_vec,rtnorm_vec
 use truncated_normal_linear, only: cholperm,cholperm_result,dmvnorm,dmvt
 use truncated_normal_core, only: prob_result, tregress_result, mvncdf, mvnqmc, mvtcdf, mvtqmc, &
  mvrandn, mvrandt, tregress
 implicit none
 private
 public :: dp, prob_result, tregress_result, cholperm_result, lnNpr, phinv, norminvp, trandn, cholperm
 public :: dmvnorm, dmvt, mvncdf, mvnqmc, mvtcdf, mvtqmc, mvrandn, mvrandt, tregress
 public::pmvnorm,pmvt,dtmvnorm,dtmvt,ptmvnorm,ptmvt,rtmvnorm,rtmvt
 public::qtnorm_vec,rtnorm_vec
contains
function pmvnorm(mu,sigma,lb,ub,B,qmc) result(r)
 real(dp),intent(in)::mu(:),sigma(:,:),lb(:),ub(:)
 integer,intent(in),optional::B
 logical,intent(in),optional::qmc
 type(prob_result)::r
 integer::n
 logical::q
 n=10000
 if(present(B))n=B
 q=.false.
 if(present(qmc))q=qmc
 if(q)then
 r=mvnqmc(lb-mu,ub-mu,sigma,n)
 else
 r=mvncdf(lb-mu,ub-mu,sigma,n)
 end if
end function
function pmvt(mu,sigma,df,lb,ub,B,qmc) result(r)
 real(dp),intent(in)::mu(:),sigma(:,:),df,lb(:),ub(:)
 integer,intent(in),optional::B
 logical,intent(in),optional::qmc
 type(prob_result)::r
 integer::n
 logical::q
 n=10000
 if(present(B))n=B
 q=.false.
 if(present(qmc))q=qmc
 if(q)then
 r=mvtqmc(lb-mu,ub-mu,sigma,df,n)
 else
 r=mvtcdf(lb-mu,ub-mu,sigma,df,n)
 end if
end function
function dtmvnorm(x,mu,sigma,lb,ub,logd,B,qmc) result(v)
 real(dp),intent(in)::x(:,:),mu(:),sigma(:,:),lb(:),ub(:)
 logical,intent(in),optional::logd,qmc
 integer,intent(in),optional::B
 real(dp)::v(size(x,1))
 type(prob_result)::r
 real(dp)::ld(size(x,1))
 integer::i
 logical::lg
 real(dp) :: pinf, ninf
 lg=.false.
 if(present(logd))lg=logd
 ld=dmvnorm(x,mu,sigma,.true.)
 pinf=ieee_value(1.0_dp,ieee_positive_inf)
 ninf=ieee_value(1.0_dp,ieee_negative_inf)
 if (.not. (all(lb==ninf) .and. all(ub==pinf))) then
  r=pmvnorm(mu,sigma,lb,ub,B,qmc)
  ld=ld-log(r%prob)
 end if
 do i=1,size(x,1)
 if(any(x(i,:)<lb).or.any(x(i,:)>ub))ld(i)=-huge(1.0_dp)
 end do
 if(lg)then
 v=ld
 else
 v=exp(ld)
 end if
end function
function dtmvt(x,mu,sigma,df,lb,ub,logd,B,qmc) result(v)
 real(dp),intent(in)::x(:,:),mu(:),sigma(:,:),df,lb(:),ub(:)
 logical,intent(in),optional::logd,qmc
 integer,intent(in),optional::B
 real(dp)::v(size(x,1))
 type(prob_result)::r
 real(dp)::ld(size(x,1))
 integer::i
 logical::lg
 real(dp) :: pinf, ninf
 pinf = ieee_value(1.0_dp, ieee_positive_inf)
 if (df == 0.0_dp .or. df == pinf) then
  v = dtmvnorm(x, mu, sigma, lb, ub, logd, B, qmc)
  return
 end if
 lg=.false.
 if(present(logd))lg=logd
 ld=dmvt(x,mu,sigma,df,.true.)
 ninf=ieee_value(1.0_dp,ieee_negative_inf)
 if (.not. (all(lb==ninf) .and. all(ub==pinf))) then
  r=pmvt(mu,sigma,df,lb,ub,B,qmc)
  ld=ld-log(r%prob)
 end if
 do i=1,size(x,1)
 if(any(x(i,:)<lb).or.any(x(i,:)>ub))ld(i)=-huge(1.0_dp)
 end do
 if(lg)then
 v=ld
 else
 v=exp(ld)
 end if
end function
function ptmvnorm(qv,mu,sigma,lb,ub,logp,B,qmc) result(v)
 real(dp),intent(in)::qv(:,:),mu(:),sigma(:,:),lb(:),ub(:)
 logical,intent(in),optional::logp,qmc
 integer,intent(in),optional::B
 real(dp)::v(size(qv,1))
 type(prob_result)::den,num
 integer::i
 logical::lg
 lg=.false.
 if(present(logp))lg=logp
 den=pmvnorm(mu,sigma,lb,ub,B,qmc)
 do i=1,size(qv,1)
  if(all(qv(i,:)>=ub))then
  v(i)=1
  else if(any(qv(i,:)<=lb))then
  v(i)=0
  else
  num=pmvnorm(mu,sigma,lb,min(ub,qv(i,:)),B,qmc)
  v(i)=min(1.0_dp,num%prob/den%prob)
  end if
 end do
 if(lg) then
  where(v>0.0_dp)
   v=log(v)
  elsewhere
   v=-huge(1.0_dp)
  end where
 end if
end function
function ptmvt(qv,mu,sigma,df,lb,ub,logp,B,qmc) result(v)
 real(dp),intent(in)::qv(:,:),mu(:),sigma(:,:),df,lb(:),ub(:)
 logical,intent(in),optional::logp,qmc
 integer,intent(in),optional::B
 real(dp)::v(size(qv,1))
 type(prob_result)::den,num
 integer::i
 logical::lg
 if (df == 0.0_dp .or. df > 350.0_dp) then
  v = ptmvnorm(qv, mu, sigma, lb, ub, logp, B, qmc)
  return
 end if
 lg=.false.
 if(present(logp))lg=logp
 den=pmvt(mu,sigma,df,lb,ub,B,qmc)
 do i=1,size(qv,1)
  if(all(qv(i,:)>=ub))then
  v(i)=1
  else if(any(qv(i,:)<=lb))then
  v(i)=0
  else
  num=pmvt(mu,sigma,df,lb,min(ub,qv(i,:)),B,qmc)
  v(i)=min(1.0_dp,num%prob/den%prob)
  end if
 end do
 if(lg) then
  where(v>0.0_dp)
   v=log(v)
  elsewhere
   v=-huge(1.0_dp)
  end where
 end if
end function
function rtmvnorm(n,mu,sigma,lb,ub) result(x)
 integer,intent(in)::n
 real(dp),intent(in)::mu(:),sigma(:,:),lb(:),ub(:)
 real(dp),allocatable::x(:,:)
 integer::d,nf,ng,i
 integer,allocatable::free_idx(:),fixed_idx(:),all_idx(:)
 real(dp),allocatable::sff(:,:),sfg(:,:),sgg(:,:),rhs(:,:),solm(:,:)
 real(dp),allocatable::delta(:),solv(:),mup(:),sigmap(:,:),xf(:,:)
 logical,allocatable::free_mask(:)

 d=size(mu)
 if(any(ub<lb)) error stop 'rtmvnorm: upper bounds must not be below lower bounds'
 allocate(free_mask(d),all_idx(d))
 all_idx=[(i,i=1,d)]
 free_mask=(ub-lb)>=1.0e-10_dp
 nf=count(free_mask)
 ng=d-nf
 if(ng==0)then
  x=mvrandn(lb,ub,sigma,n,mu)
  return
 end if
 allocate(x(n,d))
 if(nf==0)then
  do i=1,n
   x(i,:)=lb
  end do
  return
 end if
 free_idx=pack(all_idx,free_mask)
 fixed_idx=pack(all_idx,.not.free_mask)
 allocate(sff(nf,nf),sfg(nf,ng),sgg(ng,ng),delta(ng),solv(ng),mup(nf),sigmap(nf,nf))
 sff=sigma(free_idx,free_idx)
 sfg=sigma(free_idx,fixed_idx)
 sgg=sigma(fixed_idx,fixed_idx)
 delta=lb(fixed_idx)-mu(fixed_idx)
 solv=solve_real(sgg,delta)
 mup=mu(free_idx)+matmul(sfg,solv)
 allocate(rhs(ng,nf),solm(ng,nf))
 rhs=transpose(sfg)
 solm=solve_real(sgg,rhs)
 sigmap=sff-matmul(sfg,solm)
 xf=mvrandn(lb(free_idx),ub(free_idx),sigmap,n,mup)
 do i=1,n
  x(i,free_idx)=xf(i,:)
  x(i,fixed_idx)=lb(fixed_idx)
 end do
end function rtmvnorm

function rtmvt(n,mu,sigma,df,lb,ub) result(x)
 integer,intent(in)::n
 real(dp),intent(in)::mu(:),sigma(:,:),df,lb(:),ub(:)
 real(dp),allocatable::x(:,:)
 integer::d,nf,ng,i
 integer,allocatable::free_idx(:),fixed_idx(:),all_idx(:)
 real(dp),allocatable::sff(:,:),sfg(:,:),sgg(:,:),rhs(:,:),solm(:,:)
 real(dp),allocatable::delta(:),solv(:),mup(:),sigmap(:,:),xf(:,:)
 real(dp)::quad,pinf,dfp
 logical,allocatable::free_mask(:)

 pinf=ieee_value(1.0_dp,ieee_positive_inf)
 if(df==0.0_dp .or. df==pinf)then
  x=rtmvnorm(n,mu,sigma,lb,ub)
  return
 end if
 d=size(mu)
 if(any(ub<lb)) error stop 'rtmvt: upper bounds must not be below lower bounds'
 allocate(free_mask(d),all_idx(d))
 all_idx=[(i,i=1,d)]
 free_mask=(ub-lb)>=1.0e-10_dp
 nf=count(free_mask)
 ng=d-nf
 if(ng==0)then
  x=mvrandt(lb,ub,sigma,df,n,mu)
  return
 end if
 allocate(x(n,d))
 if(nf==0)then
  do i=1,n
   x(i,:)=lb
  end do
  return
 end if
 free_idx=pack(all_idx,free_mask)
 fixed_idx=pack(all_idx,.not.free_mask)
 allocate(sff(nf,nf),sfg(nf,ng),sgg(ng,ng),delta(ng),solv(ng),mup(nf),sigmap(nf,nf))
 sff=sigma(free_idx,free_idx)
 sfg=sigma(free_idx,fixed_idx)
 sgg=sigma(fixed_idx,fixed_idx)
 delta=lb(fixed_idx)-mu(fixed_idx)
 solv=solve_real(sgg,delta)
 quad=dot_product(delta,solv)
 mup=mu(free_idx)+matmul(sfg,solv)
 allocate(rhs(ng,nf),solm(ng,nf))
 rhs=transpose(sfg)
 solm=solve_real(sgg,rhs)
 dfp=df+real(ng,dp)
 sigmap=(df+quad)/dfp*(sff-matmul(sfg,solm))
 xf=mvrandt(lb(free_idx),ub(free_idx),sigmap,dfp,n,mup)
 do i=1,n
  x(i,free_idx)=xf(i,:)
  x(i,fixed_idx)=lb(fixed_idx)
 end do
end function rtmvt
end module truncated_normal_api
