module ld_importance
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface
use ld_random, only: rand_uniform, rand_mvn, rand_mv_t, rand_categorical
use ld_distributions, only: dmvn, dmvt
use ld_numerics, only: logsumexp
use ld_linalg, only: weighted_mean, weighted_covariance, make_positive_definite
implicit none
private
public :: importance_result_t, importance_normal, sir_normal, rejection_normal
public :: bayesian_bootstrap_weights, pmc_normal

type :: importance_result_t
   real(dp), allocatable :: theta(:,:), log_target(:), log_proposal(:), weights(:)
   real(dp) :: ess=0.0_dp
end type
contains
subroutine importance_normal(f,mu,sigma,n,res)
   procedure(log_target_iface)::f
   real(dp),intent(in)::mu(:),sigma(:,:)
   integer,intent(in)::n
   type(importance_result_t),intent(out)::res
   integer::i,p,info
   real(dp)::lse
   p=size(mu); allocate(res%theta(n,p),res%log_target(n),res%log_proposal(n),res%weights(n))
   do i=1,n
      call rand_mvn(mu,sigma,res%theta(i,:),info)
      res%log_target(i)=f(res%theta(i,:))
      res%log_proposal(i)=dmvn(res%theta(i,:),mu,sigma,.true.)
   end do
   lse=logsumexp(res%log_target-res%log_proposal); res%weights=exp(res%log_target-res%log_proposal-lse)
   res%ess=1.0_dp/sum(res%weights**2)
end subroutine importance_normal

subroutine sir_normal(f,mu,sigma,n,theta_out)
   procedure(log_target_iface)::f
   real(dp),intent(in)::mu(:),sigma(:,:)
   integer,intent(in)::n
   real(dp),intent(out)::theta_out(:,:)
   type(importance_result_t)::res
   integer::i,k
   call importance_normal(f,mu,sigma,n,res)
   do i=1,size(theta_out,1); k=rand_categorical(res%weights); theta_out(i,:)=res%theta(k,:); end do
end subroutine sir_normal

subroutine rejection_normal(f,mu,sigma,logc,n_try,accepted,nacc)
   procedure(log_target_iface)::f
   real(dp),intent(in)::mu(:),sigma(:,:),logc
   integer,intent(in)::n_try
   real(dp),intent(out)::accepted(:,:)
   integer,intent(out)::nacc
   integer::i,info
   real(dp)::x(size(mu)),lp,lg
   nacc=0
   do i=1,n_try
      call rand_mvn(mu,sigma,x,info); lp=f(x); lg=dmvn(x,mu,sigma,.true.)
      if(log(rand_uniform())<min(0.0_dp,lp-lg-logc)) then
         if(nacc<size(accepted,1)) then; nacc=nacc+1; accepted(nacc,:)=x; end if
      end if
   end do
end subroutine rejection_normal

subroutine bayesian_bootstrap_weights(nobs,nsim,w)
   integer,intent(in)::nobs,nsim
   real(dp),intent(out)::w(nobs,nsim)
   integer::s,i,j
   real(dp)::u(max(1,nobs-1)),tmp
   if(nobs==1) then; w=1.0_dp; return; end if
   do s=1,nsim
      do i=1,nobs-1; u(i)=rand_uniform(); end do
      do i=1,nobs-2; do j=i+1,nobs-1; if(u(j)<u(i)) then; tmp=u(i);u(i)=u(j);u(j)=tmp; end if; end do; end do
      w(1,s)=u(1); do i=2,nobs-1; w(i,s)=u(i)-u(i-1); end do; w(nobs,s)=1.0_dp-u(nobs-1)
   end do
end subroutine bayesian_bootstrap_weights

subroutine pmc_normal(f,mu0,cov0,n_particles,n_iter,mu,cov,ess_history)
   procedure(log_target_iface)::f
   real(dp),intent(in)::mu0(:),cov0(:,:)
   integer,intent(in)::n_particles,n_iter
   real(dp),intent(out)::mu(:),cov(:,:),ess_history(:)
   type(importance_result_t)::res
   integer::it
   mu=mu0; cov=cov0
   do it=1,n_iter
      call importance_normal(f,mu,cov,n_particles,res)
      mu=weighted_mean(res%theta,res%weights)
      cov=weighted_covariance(res%theta,res%weights,mu)
      cov=cov+1e-10_dp*identity(size(mu)); call make_positive_definite(cov); if(it<=size(ess_history)) ess_history(it)=res%ess
   end do
contains
   pure function identity(n) result(a)
      integer,intent(in)::n; real(dp)::a(n,n); integer::i
      a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
   end function identity
end subroutine pmc_normal
end module ld_importance
