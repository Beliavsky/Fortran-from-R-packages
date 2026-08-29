! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_gel
use r_compat, only: dp
use gmm_linalg, only: solve_linear, least_squares, quadratic_inverse
implicit none
private
public :: gel_rho, gel_lambda, gel_implied_prob, gel_objective, gel_result_t

integer, parameter, public :: GEL_EL=1, GEL_ET=2, GEL_CUE=3, GEL_ETEL=4, GEL_HD=5, GEL_ETHD=6, GEL_RCUE=7

type :: gel_result_t
   real(dp), allocatable :: lambda(:), prob(:)
   real(dp) :: objective=0.0_dp
   integer :: convergence=0, iterations=0
end type gel_result_t

contains

pure function gel_rho(gt,lambda,derive,type,k) result(r)
real(dp),intent(in)::gt(:,:),lambda(:)
integer,intent(in)::derive,type
real(dp),intent(in),optional::k
real(dp)::r(size(gt,1)),gml(size(gt,1)),kk,w(size(gt,1)),sw
kk=1.0_dp
if(present(k)) kk=k
gml=matmul(gt,lambda)*kk
select case(derive)
case(0)
   select case(type)
   case(GEL_EL)
      where(1.0_dp-gml>0.0_dp)
      r=log(1.0_dp-gml)
      elsewhere
      r=-huge(1.0_dp)
      end where
   case(GEL_ET)
      r=-exp(gml)
   case(GEL_CUE,GEL_RCUE)
      r=-gml-0.5_dp*gml*gml
   case(GEL_HD)
      r=-2.0_dp/(1.0_dp-gml/2.0_dp)
   case(GEL_ETEL)
      w=exp(gml)
      sw=sum(w)
      if(sw<=0) then
      r=huge(1.0_dp)
      else
      w=w/sw
      r=-log(max(w*real(size(gt,1),dp),tiny(1.0_dp)))
      end if
   case(GEL_ETHD)
      w=exp(gml)
      sw=sum(w)
      if(sw<=0) then
      r=huge(1.0_dp)
      else
      w=w/sw
      r=(sqrt(w)-1.0_dp/sqrt(real(size(gt,1),dp)))**2
      end if
   end select
case(1)
   select case(type)
   case(GEL_EL); r=-1.0_dp/(1.0_dp-gml)
   case(GEL_ET,GEL_ETEL,GEL_ETHD); r=-exp(gml)
   case(GEL_CUE,GEL_RCUE); r=-1.0_dp-gml
   case(GEL_HD); r=-1.0_dp/(1.0_dp-gml/2.0_dp)**2
   end select
case(2)
   select case(type)
   case(GEL_EL); r=-1.0_dp/(1.0_dp-gml)**2
   case(GEL_ET,GEL_ETEL,GEL_ETHD); r=-exp(gml)
   case(GEL_CUE,GEL_RCUE); r=-1.0_dp
   case(GEL_HD); r=-1.0_dp/(1.0_dp-gml/2.0_dp)**3
   end select
end select
end function gel_rho

pure function rho0_at_zero(type) result(v)
integer,intent(in)::type
real(dp)::v
select case(type)
case(GEL_EL); v=0.0_dp
case(GEL_ET); v=-1.0_dp
case(GEL_CUE,GEL_RCUE); v=0.0_dp
case(GEL_HD); v=-2.0_dp
case(GEL_ETEL,GEL_ETHD); v=0.0_dp
end select
end function rho0_at_zero

function gel_objective(gt,lambda,type,k) result(obj)
real(dp),intent(in)::gt(:,:),lambda(:)
integer,intent(in)::type
real(dp),intent(in),optional::k
real(dp)::obj,kk
real(dp)::r(size(gt,1))
kk=1.0_dp
if(present(k)) kk=k
r=gel_rho(gt,lambda,0,type,kk)
if(type==GEL_ETHD) then
   obj=sum(r)
else
   obj=sum(r)/real(size(r),dp)-rho0_at_zero(type)
end if
end function gel_objective

subroutine gel_lambda(gt,type,res,tol_lambda,maxit,tol_obj,k)
real(dp),intent(in)::gt(:,:)
integer,intent(in)::type
class(gel_result_t),intent(out)::res
real(dp),intent(in),optional::tol_lambda,tol_obj,k
integer,intent(in),optional::maxit
real(dp)::tol,tobj,kk
integer::mi,solve_type
allocate(res%lambda(size(gt,2)),res%prob(size(gt,1)))
res%lambda=0
res%prob=1.0_dp/size(gt,1)
tol=1.0e-8_dp
if(present(tol_lambda)) tol=tol_lambda
tobj=1.0e-8_dp
if(present(tol_obj)) tobj=tol_obj
mi=100
if(present(maxit)) mi=maxit
kk=1.0_dp
if(present(k)) kk=k
select case(type)
case(GEL_EL)
   call lambda_wu(gt,res%lambda,res%convergence,res%iterations,tol,mi)
case(GEL_CUE)
   call lambda_cue(gt,res%lambda,res%convergence)
case(GEL_RCUE)
   call lambda_rcue(gt,res%lambda,res%prob,res%convergence,res%iterations,mi,kk)
case default
   solve_type=type
   if(type==GEL_ETEL .or. type==GEL_ETHD) solve_type=GEL_ET
   call lambda_newton(gt,solve_type,res%lambda,res%convergence,res%iterations,tol,tobj,mi,kk)
end select
if(type/=GEL_RCUE) res%prob=gel_implied_prob(gt,res%lambda,type,kk)
res%objective=gel_objective(gt,res%lambda,type,kk)
end subroutine gel_lambda

subroutine lambda_wu(gt,lambda,conv,iters,tol,maxit)
real(dp),intent(in)::gt(:,:),tol
real(dp),intent(out)::lambda(:)
integer,intent(out)::conv,iters
integer,intent(in)::maxit
real(dp)::m(size(gt,2)),d1(size(gt,2)),tmp(size(gt,1)),den(size(gt,1))
real(dp)::dd(size(gt,2),size(gt,2)),dif
real(dp),allocatable::step(:)
integer::i,info
m=0.0_dp
conv=1
DO i=1,maxit
   den=1.0_dp+matmul(gt,m)
   if(any(den<=0.0_dp)) exit
   d1=matmul(transpose(gt),1.0_dp/den)
   dd=-matmul(transpose(gt),gt*spread(1.0_dp/den**2,2,size(gt,2)))
   call solve_linear(dd,d1,step,info)
   if(info/=0) exit
   dif=maxval(abs(step))
   do
      tmp=1.0_dp+matmul(gt,m-step)
      if(minval(tmp)>0.0_dp) exit
      step=step/2.0_dp
      if(maxval(abs(step))<epsilon(1.0_dp)) exit
   end do
   m=m-step
   if(dif<=tol) then
   conv=0
   exit
   end if
END DO
iters=i
if(conv/=0) m=0.0_dp
lambda=-m
end subroutine lambda_wu

subroutine lambda_newton(gt,type,lambda,conv,iters,tol,tobj,maxit,k)
real(dp),intent(in)::gt(:,:),tol,tobj,k
integer,intent(in)::type,maxit
real(dp),intent(out)::lambda(:)
integer,intent(out)::conv,iters
real(dp)::r1(size(gt,1)),r2(size(gt,1)),f(size(gt,2)),jmat(size(gt,2),size(gt,2))
real(dp),allocatable::step(:)
integer::i,info
lambda=0.0_dp
conv=1
do i=1,maxit
   r1=gel_rho(gt,lambda,1,type,k)
   r2=gel_rho(gt,lambda,2,type,k)
   f=-sum(gt*spread(r1,2,size(gt,2)),dim=1)/real(size(gt,1),dp)
   if(sum(abs(f))<tobj) then
   conv=0
   exit
   end if
   jmat=matmul(transpose(gt*spread(r2,2,size(gt,2))),gt)
   call solve_linear(jmat,f,step,info)
   if(info/=0) exit
   if(sum(abs(step))<tol) then
   conv=0
   exit
   end if
   lambda=lambda+step
end do
iters=i
end subroutine lambda_newton

subroutine lambda_cue(gt,lambda,conv)
real(dp),intent(in)::gt(:,:)
real(dp),intent(out)::lambda(:)
integer,intent(out)::conv
real(dp)::minus_one(size(gt,1))
real(dp),allocatable::lamtmp(:)
integer::info
minus_one=-1.0_dp
call least_squares(gt,minus_one,lamtmp,info)
if(info==0) then
lambda=lamtmp
else
lambda=0.0_dp
end if
conv=merge(0,1,info==0)
end subroutine lambda_cue

subroutine lambda_rcue(gt,lambda,prob,conv,iters,maxit,k)
real(dp),intent(in)::gt(:,:),k
real(dp),intent(out)::lambda(:),prob(:)
integer,intent(out)::conv,iters
integer,intent(in)::maxit
logical::keep(size(gt,1)),newkeep(size(gt,1))
integer::i,n1,info
real(dp),allocatable::sub(:,:),lam(:)
real(dp)::raw(size(gt,1))
keep=.true.
conv=0
lambda=0.0_dp
prob=1.0_dp/size(gt,1)
do i=1,maxit
   n1=count(keep)
   if(n1<size(gt,2)+1) then
   conv=2
   exit
   end if
   allocate(sub(n1,size(gt,2)))
   sub=pack_rows(gt,keep)
   allocate(lam(size(gt,2)))
   call lambda_cue(sub,lam,info)
   deallocate(sub)
   if(info/=0) then
   conv=3
   deallocate(lam)
   exit
   end if
   lambda=lam
   deallocate(lam)
   raw=1.0_dp+k*matmul(gt,lambda)
   newkeep=raw>0.0_dp
   if(all(newkeep.eqv.keep)) exit
   keep=newkeep
end do
iters=i
if(i>maxit) conv=1
if(conv/=0) then
lambda=0
prob=1.0_dp/size(gt,1)
else
prob=max(raw,0.0_dp)
prob=prob/sum(prob)
end if
end subroutine lambda_rcue

pure function pack_rows(a,mask) result(b)
real(dp),intent(in)::a(:,:)
logical,intent(in)::mask(:)
real(dp)::b(count(mask),size(a,2))
integer::i,j
j=0
do i=1,size(a,1)
   if(mask(i)) then
   j=j+1
   b(j,:)=a(i,:)
   end if
end do
end function pack_rows

pure function gel_implied_prob(gt,lambda,type,k) result(pt)
real(dp),intent(in)::gt(:,:),lambda(:),k
integer,intent(in)::type
real(dp)::pt(size(gt,1)),eps
integer::t
t=type
if(t==GEL_ETEL .or. t==GEL_ETHD) t=GEL_ET
pt=-gel_rho(gt,lambda,1,t,k)/real(size(gt,1),dp)
if(t==GEL_CUE) then
   eps=-real(size(pt),dp)*min(minval(pt),0.0_dp)
   pt=(pt+eps/real(size(pt),dp))/(1.0_dp+eps)
else if(t==GEL_RCUE) then
   where(pt<0.0_dp) pt=0.0_dp
end if
if(sum(pt)>0.0_dp) pt=pt/sum(pt)
end function gel_implied_prob

end module gmm_gel
