! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_ategel_fit
use r_compat, only: dp
use gmm_ate, only: ate_moments,ate_gradient,ATE_BAL,ATE_LINEAR
use gmm_gel, only: GEL_EL
use gmm_gel_fit, only: gel_fit_result_t,gel_fit
implicit none
private
public :: ategel_fit
integer,save :: s_k=0,s_mom_type=ATE_BAL,s_family=ATE_LINEAR
real(dp),allocatable,save :: s_w(:,:),s_pop(:)
contains
subroutine ategel_fit(x,k,theta0,res,mom_type,family,type,w,pop_mom,tol,maxit,maxiterlam)
real(dp),intent(in)::x(:,:),theta0(:)
integer,intent(in)::k
class(gel_fit_result_t),intent(out)::res
integer,intent(in),optional::mom_type,family,type,maxit,maxiterlam
real(dp),intent(in),optional::w(:,:),pop_mom(:),tol
integer::mt,fam,tp
s_k=k
mt=ATE_BAL
if(present(mom_type))mt=mom_type
s_mom_type=mt
fam=ATE_LINEAR
if(present(family))fam=family
s_family=fam
tp=GEL_EL
if(present(type))tp=type
if(allocated(s_w))deallocate(s_w)
if(allocated(s_pop))deallocate(s_pop)
if(present(w))s_w=w
if(present(pop_mom))s_pop=pop_mom
if(present(tol))then
 if(present(maxit).and.present(maxiterlam))then
  call gel_fit(ate_callback,x,theta0,tp,res,gradient=ate_grad_callback,tol=tol,maxit=maxit,maxiterlam=maxiterlam)
 else if(present(maxit))then
  call gel_fit(ate_callback,x,theta0,tp,res,gradient=ate_grad_callback,tol=tol,maxit=maxit)
 else
  call gel_fit(ate_callback,x,theta0,tp,res,gradient=ate_grad_callback,tol=tol)
 end if
else
 call gel_fit(ate_callback,x,theta0,tp,res,gradient=ate_grad_callback)
end if
if(allocated(s_w))deallocate(s_w)
if(allocated(s_pop))deallocate(s_pop)
end subroutine ategel_fit

pure function ate_callback(theta,data) result(gt)
real(dp),intent(in)::theta(:),data(:,:)
real(dp),allocatable::gt(:,:)
if(allocated(s_w).and.allocated(s_pop))then
 call ate_moments(theta,data,s_k,s_mom_type,s_family,gt,w=s_w,pop_mom=s_pop)
else if(allocated(s_w))then
 call ate_moments(theta,data,s_k,s_mom_type,s_family,gt,w=s_w)
else if(allocated(s_pop))then
 call ate_moments(theta,data,s_k,s_mom_type,s_family,gt,pop_mom=s_pop)
else
 call ate_moments(theta,data,s_k,s_mom_type,s_family,gt)
end if
end function ate_callback

pure function ate_grad_callback(theta,data) result(g)
real(dp),intent(in)::theta(:),data(:,:)
real(dp),allocatable::g(:,:)
if(allocated(s_w).and.allocated(s_pop))then
 call ate_gradient(theta,data,s_k,s_mom_type,s_family,g,w=s_w,pop_mom=s_pop)
else if(allocated(s_w))then
 call ate_gradient(theta,data,s_k,s_mom_type,s_family,g,w=s_w)
else if(allocated(s_pop))then
 call ate_gradient(theta,data,s_k,s_mom_type,s_family,g,pop_mom=s_pop)
else
 call ate_gradient(theta,data,s_k,s_mom_type,s_family,g)
end if
end function ate_grad_callback
end module gmm_ategel_fit
