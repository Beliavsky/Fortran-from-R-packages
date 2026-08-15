module ld_mcmc_completion
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface
use ld_random, only: rand_uniform, rand_normal, rand_mvn
use ld_numerics, only: numerical_gradient
use ld_linalg, only: sample_covariance, make_positive_definite, chol_lower
use ld_mcmc, only: mcmc_result_t, mwg_sample
use ld_mcmc_advanced, only: dram_sample
implicit none
private
public :: admg_sample, adaptive_hmc_sample, drm_sample, gibbs_sample
public :: inca_sample, refractive_sample, smwg_sample, samwg_sample
public :: usmwg_sample, usamwg_sample, reversible_jump_selection_sample
public :: rj_selection_result_t, gibbs_draw_iface

abstract interface
   subroutine gibbs_draw_iface(index,state)
      import dp
      integer, intent(in) :: index
      real(dp), intent(inout) :: state(:)
   end subroutine gibbs_draw_iface
end interface

type :: rj_selection_result_t
   real(dp), allocatable :: chain(:,:), logp(:)
   logical, allocatable :: active(:,:)
   real(dp) :: acceptance_rate=0.0_dp
   integer :: accepted=0, proposed=0
end type rj_selection_result_t

contains

subroutine init_result(res,nkeep,p)
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in) :: nkeep,p
   allocate(res%chain(nkeep,p),res%logp(nkeep))
   res%chain=0.0_dp; res%logp=0.0_dp
   res%acceptance_rate=0.0_dp; res%accepted=0; res%proposed=0
end subroutine init_result

subroutine random_order(ord)
   integer, intent(out) :: ord(:)
   integer :: i,j,t
   do i=1,size(ord); ord(i)=i; end do
   do i=size(ord),2,-1
      j=1+int(rand_uniform()*real(i,dp))
      t=ord(i); ord(i)=ord(j); ord(j)=t
   end do
end subroutine random_order

subroutine admg_sample(f,x0,initial_cov,n_iter,burn,thin,res,periodicity,n_prior,final_cov)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),initial_cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: periodicity,n_prior
   real(dp), allocatable, intent(out), optional :: final_cov(:,:)
   integer :: p,it,j,k,nkeep,period,n0,ord(size(x0)),info
   real(dp) :: x(size(x0)),prop(size(x0)),lp,lpp,sd(size(x0)),acc(size(x0))
   real(dp), allocatable :: history(:,:),vc(:,:)
   p=size(x0); period=50; if(present(periodicity)) period=max(1,periodicity)
   n0=10; if(present(n_prior)) n0=max(1,n_prior)
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result(res,nkeep,p)
   allocate(history(n_iter+1,p),vc(p,p)); vc=initial_cov; call make_positive_definite(vc)
   sd=sqrt(max(diagonal(vc),1.0e-10_dp)); x=x0; lp=f(x); history(1,:)=x; acc=0.0_dp; k=0
   do it=1,n_iter
      if(it>10) sd=sqrt(max(diagonal(vc),1.0e-10_dp))*exp(diagonal(vc)*(acc/real(it-1,dp)-0.30_dp))
      call random_order(ord)
      do j=1,p
         prop=x; prop(ord(j))=x(ord(j))+sqrt(0.01_dp+sd(ord(j))**2)*rand_normal(); lpp=f(prop)
         res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then
            x=prop; lp=lpp; res%accepted=res%accepted+1; acc(ord(j))=acc(ord(j))+1.0_dp
         end if
      end do
      history(it+1,:)=x
      if(mod(it,period)==0 .and. it+n0>p+1) then
         call sample_covariance(history(1:it+1,:),vc)
         vc=vc+1.0e-5_dp*identity(p); call make_positive_definite(vc)
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_cov)) then; allocate(final_cov(p,p)); final_cov=vc; end if
end subroutine admg_sample

subroutine adaptive_hmc_sample(f,x0,epsilon,n_leapfrog,n_iter,burn,thin,res,periodicity,final_epsilon)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),epsilon(:)
   integer, intent(in) :: n_leapfrog,n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: periodicity
   real(dp), allocatable, intent(out), optional :: final_epsilon(:)
   integer :: p,it,l,j,k,nkeep,period,winacc
   real(dp) :: x(size(x0)),q(size(x0)),mom(size(x0)),mom0(size(x0))
   real(dp) :: g(size(x0)),eps(size(x0)),lp,lpp,h0,h1
   p=size(x0); period=10; if(present(periodicity)) period=max(1,periodicity)
   eps=epsilon; nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result(res,nkeep,p)
   x=x0; lp=f(x); k=0; winacc=0
   do it=1,n_iter
      do j=1,p; mom0(j)=rand_normal(); end do
      mom=mom0; q=x; call numerical_gradient(f,q,g); mom=mom+0.5_dp*eps*g
      do l=1,max(1,n_leapfrog)
         q=q+eps*mom; call numerical_gradient(f,q,g)
         if(l<n_leapfrog) mom=mom+eps*g
      end do
      mom=mom+0.5_dp*eps*g; lpp=f(q)
      h0=-lp+0.5_dp*dot_product(mom0,mom0); h1=-lpp+0.5_dp*dot_product(mom,mom)
      res%proposed=res%proposed+1
      if(log(rand_uniform())<min(0.0_dp,h0-h1)) then
         x=q; lp=lpp; res%accepted=res%accepted+1; winacc=winacc+1
      end if
      if(mod(it,period)==0 .and. it<=burn) then
         if(real(winacc,dp)/real(period,dp)<=0.1_dp) eps=0.8_dp*eps
         if(real(winacc,dp)/real(period,dp)>0.7_dp) eps=1.2_dp*eps
         winacc=0
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_epsilon)) then; allocate(final_epsilon(p)); final_epsilon=eps; end if
end subroutine adaptive_hmc_sample

subroutine drm_sample(f,x0,cov,n_iter,burn,thin,res,scale_factor)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   real(dp), intent(in), optional :: scale_factor
   real(dp) :: sf
   sf=0.5_dp; if(present(scale_factor)) sf=scale_factor
   call dram_sample(f,x0,cov,n_iter,burn,thin,res,adapt_start=n_iter+1,periodicity=max(1,n_iter),scale_factor=sf)
end subroutine drm_sample

subroutine gibbs_sample(f,draw_conditional,x0,n_iter,burn,thin,res,random_scan)
   procedure(log_target_iface) :: f
   procedure(gibbs_draw_iface) :: draw_conditional
   real(dp), intent(in) :: x0(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   logical, intent(in), optional :: random_scan
   integer :: p,it,j,k,nkeep,ord(size(x0))
   logical :: rs
   real(dp) :: x(size(x0))
   p=size(x0); rs=.true.; if(present(random_scan)) rs=random_scan
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result(res,nkeep,p); x=x0; k=0
   do it=1,n_iter
      if(rs) then; call random_order(ord); else; do j=1,p; ord(j)=j; end do; end if
      do j=1,p; call draw_conditional(ord(j),x); end do
      res%proposed=res%proposed+p; res%accepted=res%accepted+p
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=f(x)
      end if
   end do
   res%acceptance_rate=1.0_dp
end subroutine gibbs_sample

subroutine inca_sample(f,chains,initial_cov,n_iter,burn,thin,out,acceptance,periodicity,adaptive_start)
   procedure(log_target_iface) :: f
   real(dp), intent(inout) :: chains(:,:)
   real(dp), intent(in) :: initial_cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   real(dp), allocatable, intent(out) :: out(:,:,:)
   real(dp), intent(out) :: acceptance
   integer, intent(in), optional :: periodicity,adaptive_start
   integer :: nc,p,it,i,k,nkeep,period,astart,info,acc,propn
   real(dp) :: cov(size(chains,2),size(chains,2)),step(size(chains,2)),prop(size(chains,2))
   real(dp) :: lp(size(chains,1)),lpp,scale
   real(dp), allocatable :: pool(:,:)
   nc=size(chains,1); p=size(chains,2); period=50; if(present(periodicity)) period=max(2,periodicity)
   astart=100; if(present(adaptive_start)) astart=max(2,adaptive_start)
   nkeep=max(0,(n_iter-burn)/max(1,thin)); allocate(out(nkeep,nc,p)); out=0.0_dp
   cov=initial_cov; call make_positive_definite(cov); scale=2.38_dp**2/real(max(1,p),dp)
   do i=1,nc; lp(i)=f(chains(i,:)); end do
   k=0; acc=0; propn=0
   do it=1,n_iter
      do i=1,nc
         call rand_mvn(0.0_dp*step,cov,step,info); if(info/=0) cycle
         prop=chains(i,:)+step; lpp=f(prop); propn=propn+1
         if(log(rand_uniform())<lpp-lp(i)) then; chains(i,:)=prop; lp(i)=lpp; acc=acc+1; end if
      end do
      if(it>=astart .and. mod(it,period)==0) then
         allocate(pool(nc,p)); pool=chains; call sample_covariance(pool,cov); deallocate(pool)
         cov=scale*cov+1.0e-9_dp*identity(p); call make_positive_definite(cov)
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then; k=k+1; out(k,:,:)=chains; end if
   end do
   acceptance=real(acc,dp)/real(max(1,propn),dp)
end subroutine inca_sample

subroutine refractive_sample(f,x0,width,refractive_index,n_steps,n_iter,burn,thin,res,adapt)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),width,refractive_index
   integer, intent(in) :: n_steps,n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   logical, intent(in), optional :: adapt
   integer :: p,it,i,j,k,nkeep
   real(dp) :: x(size(x0)),prop(size(x0)),mom(size(x0)),g(size(x0)),u(size(x0))
   real(dp) :: lp,lpp,w,ridx,r1,r2,c1,c2,nm,ng,jac,loga,astar
   logical :: ad
   p=size(x0); w=width; ridx=max(refractive_index,1.0e-8_dp); astar=0.65_dp
   ad=.false.; if(present(adapt)) ad=adapt
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result(res,nkeep,p); x=x0; lp=f(x); k=0
   do it=1,n_iter
      prop=x; do j=1,p; mom(j)=rand_normal(); end do; jac=1.0_dp
      call numerical_gradient(f,prop,g)
      do i=1,max(1,n_steps)
         ng=sqrt(max(dot_product(g,g),tiny(1.0_dp))); nm=sqrt(max(dot_product(mom,mom),tiny(1.0_dp)))
         if(dot_product(mom,g)>0.0_dp) then; u=g/ng; r1=1.0_dp; r2=ridx
         else; u=-g/ng; r1=ridx; r2=1.0_dp; end if
         c1=dot_product(mom,u)/nm; c2=1.0_dp-(r1/r2)**2*(1.0_dp-c1*c1)
         if(c2<0.0_dp) then
            mom=mom-2.0_dp*dot_product(mom,u)*u
         else
            c2=sqrt(c2)
            mom=(r1/r2)*mom-nm*((r1/r2)*c1-c2)*u
            if(abs(c2)>tiny(1.0_dp)) jac=jac*(r1/r2)**(p-1)*(c1/c2)
         end if
         prop=prop+w*mom; lpp=f(prop); call numerical_gradient(f,prop,g)
      end do
      loga=lpp-lp+log(max(abs(jac),tiny(1.0_dp))); res%proposed=res%proposed+1
      if(log(rand_uniform())<loga) then
         x=prop; lp=lpp; res%accepted=res%accepted+1
         if(ad .and. it<=burn) w=w+(w/(astar*(1.0_dp-astar)))*(1.0_dp-astar)/real(it,dp)
      else if(ad .and. it<=burn) then
         w=abs(w-(w/(astar*(1.0_dp-astar)))*astar/real(it,dp))
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine refractive_sample

subroutine smwg_sample(f,x0,proposal_sd,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),proposal_sd(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   call mwg_sample(f,x0,proposal_sd,n_iter,burn,thin,res)
end subroutine smwg_sample

subroutine samwg_sample(f,x0,proposal_sd,n_iter,burn,thin,res,periodicity,final_sd)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),proposal_sd(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: periodicity
   real(dp), allocatable, intent(out), optional :: final_sd(:)
   integer :: p,it,j,k,nkeep,period,ord(size(x0))
   real(dp) :: x(size(x0)),prop(size(x0)),sd(size(x0)),lp,lpp,acc(size(x0)),eta
   p=size(x0); period=50; if(present(periodicity)) period=max(1,periodicity)
   sd=proposal_sd; acc=0.0_dp; nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      call random_order(ord)
      do j=1,p
         prop=x; prop(ord(j))=x(ord(j))+sd(ord(j))*rand_normal(); lpp=f(prop); res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then
            x=prop; lp=lpp; res%accepted=res%accepted+1; acc(ord(j))=acc(ord(j))+1.0_dp
         end if
      end do
      if(mod(it,period)==0 .and. it<=burn) then
         eta=1.0_dp/min(100.0_dp,sqrt(real(it,dp)))
         do j=1,p
            if(acc(j)/real(it,dp)>0.44_dp) then; sd(j)=sd(j)*exp(eta); else; sd(j)=sd(j)*exp(-eta); end if
         end do
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_sd)) then; allocate(final_sd(p)); final_sd=sd; end if
end subroutine samwg_sample

subroutine usmwg_sample(f,baseline,dynamic,x0,proposal_sd,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: baseline(:,:),x0(:),proposal_sd(:)
   integer, intent(in) :: dynamic(:),n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   call updating_mwg_core(f,baseline,dynamic,x0,proposal_sd,n_iter,burn,thin,res,.false.)
end subroutine usmwg_sample

subroutine usamwg_sample(f,baseline,dynamic,x0,proposal_sd,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: baseline(:,:),x0(:),proposal_sd(:)
   integer, intent(in) :: dynamic(:),n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   call updating_mwg_core(f,baseline,dynamic,x0,proposal_sd,n_iter,burn,thin,res,.true.)
end subroutine usamwg_sample

subroutine updating_mwg_core(f,baseline,dynamic,x0,proposal_sd,n_iter,burn,thin,res,adaptive)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: baseline(:,:),x0(:),proposal_sd(:)
   integer, intent(in) :: dynamic(:),n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   logical, intent(in) :: adaptive
   integer :: p,it,j,k,nkeep,row,idx,ord(size(dynamic))
   real(dp) :: x(size(x0)),prop(size(x0)),sd(size(x0)),lp,lpp,acc(size(x0)),eta
   p=size(x0); sd=proposal_sd; acc=0.0_dp; nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result(res,nkeep,p); k=0
   x=x0
   do it=1,n_iter
      row=1+mod(it-1,size(baseline,1)); prop=baseline(row,:); prop(dynamic)=x(dynamic); x=prop; lp=f(x)
      call random_order(ord)
      do j=1,size(dynamic)
         idx=dynamic(ord(j)); prop=x; prop(idx)=x(idx)+sd(idx)*rand_normal(); lpp=f(prop); res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; acc(idx)=acc(idx)+1.0_dp; end if
      end do
      if(adaptive .and. mod(it,50)==0 .and. it<=burn) then
         eta=1.0_dp/min(100.0_dp,sqrt(real(it,dp)))
         do j=1,size(dynamic); idx=dynamic(j)
            if(acc(idx)/real(it,dp)>0.44_dp) then; sd(idx)=sd(idx)*exp(eta); else; sd(idx)=sd(idx)*exp(-eta); end if
         end do
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine updating_mwg_core

subroutine reversible_jump_selection_sample(f,x0,selectable,initial_active,max_selected,parm_prob,bin_prob, &
      n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),parm_prob(:),bin_prob
   logical, intent(in) :: selectable(:),initial_active(:)
   integer, intent(in) :: max_selected,n_iter,burn,thin
   type(rj_selection_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,j,m,actn,ord(size(x0)),cand(size(x0)),ncand
   real(dp) :: x(size(x0)),prop(size(x0)),saved(size(x0)),theta(size(x0)),lp,lpp,lambda
   real(dp) :: prior0,prior1,loga
   logical :: active(size(x0)),aprop(size(x0))
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin))
   allocate(res%chain(nkeep,p),res%active(nkeep,p),res%logp(nkeep)); res%chain=0.0_dp; res%active=.false.; res%logp=0.0_dp
   active=initial_active; where(.not.selectable) active=.true.; x=x0; where(.not.active) x=0.0_dp; saved=x; lp=f(x); k=0
   do it=1,n_iter
      ncand=0; do j=1,p; if(selectable(j)) then; ncand=ncand+1; cand(ncand)=j; end if; end do
      if(ncand==0) exit
      j=cand(1+int(rand_uniform()*real(ncand,dp))); aprop=active
      if(count(active .and. selectable)<max_selected .or. active(j)) aprop(j)=.not.active(j)
      prior0=selection_logprior(active,selectable,max_selected,parm_prob,bin_prob)
      prior1=selection_logprior(aprop,selectable,max_selected,parm_prob,bin_prob)
      do m=1,p; theta(m)=rand_normal(); end do; theta=theta/sqrt(max(dot_product(theta,theta),tiny(1.0_dp)))
      lambda=rand_uniform(); call random_order(ord)
      do m=1,p
         if(.not.active(ord(m))) cycle
         prop=x; prop(ord(m))=x(ord(m))+lambda*theta(ord(m)); lpp=f(prop); res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; saved(ord(m))=x(ord(m)); res%accepted=res%accepted+1; end if
      end do
      prop=x
      if(aprop(j)) then
         if(saved(j)==0.0_dp) saved(j)=lambda*theta(j)
         prop(j)=saved(j)+lambda*theta(j)
      else
         if(x(j)/=0.0_dp) saved(j)=x(j); prop(j)=0.0_dp
      end if
      lpp=f(prop); loga=lpp-lp+prior1-prior0; res%proposed=res%proposed+1
      if(log(rand_uniform())<loga) then; x=prop; lp=lpp; active=aprop; res%accepted=res%accepted+1; end if
      actn=count(active)
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%active(k,:)=active; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine reversible_jump_selection_sample

function selection_logprior(active,selectable,max_selected,parm_prob,bin_prob) result(v)
   logical, intent(in) :: active(:),selectable(:)
   integer, intent(in) :: max_selected
   real(dp), intent(in) :: parm_prob(:),bin_prob
   real(dp) :: v,pb,pp
   integer :: j,n,k
   v=0.0_dp; n=count(selectable); k=count(active .and. selectable); pb=min(max(bin_prob,1.0e-12_dp),1.0_dp-1.0e-12_dp)
   do j=1,size(active)
      if(.not.selectable(j)) cycle
      pp=min(max(parm_prob(j),1.0e-12_dp),1.0_dp-1.0e-12_dp)
      if(active(j)) then; v=v+log(pp); else; v=v+log(1.0_dp-pp); end if
   end do
   if(k<=max_selected) v=v+log_binom_pmf(k,max_selected,pb)
end function selection_logprior

function log_binom_pmf(k,n,p) result(v)
   integer, intent(in) :: k,n
   real(dp), intent(in) :: p
   real(dp) :: v
   if(k<0 .or. k>n) then; v=-huge(1.0_dp); return; end if
   v=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp)) &
      +real(k,dp)*log(p)+real(n-k,dp)*log(1.0_dp-p)
end function log_binom_pmf

pure function diagonal(a) result(d)
   real(dp), intent(in) :: a(:,:)
   real(dp) :: d(min(size(a,1),size(a,2)))
   integer :: i
   do i=1,size(d); d(i)=a(i,i); end do
end function diagonal

pure function identity(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
end function identity

end module ld_mcmc_completion
