module ld_mcmc_catalog
use ld_kinds, only: dp, pi
use ld_interfaces, only: log_target_iface
use ld_random, only: rand_uniform, rand_normal, rand_mvn
use ld_linalg, only: make_positive_definite, sample_covariance, chol_lower
use ld_numerics, only: numerical_gradient
use ld_mcmc, only: mcmc_result_t
implicit none
private
public :: independence_metropolis_sample, multiple_try_metropolis_sample
public :: hit_and_run_sample, charm_sample, adaptive_mwg_sample
public :: adaptive_mixture_metropolis_sample, sgld_sample, griddy_gibbs_sample
public :: adaptive_griddy_gibbs_sample, ohss_sample, uess_sample, afss_sample
public :: mcmcmc_sample, random_dive_sample, rss_sample, tempered_hmc_sample

contains

subroutine init_result_cat(res,nkeep,p)
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in) :: nkeep,p
   allocate(res%chain(nkeep,p),res%logp(nkeep))
   res%chain=0.0_dp; res%logp=0.0_dp
   res%acceptance_rate=0.0_dp; res%proposed=0; res%accepted=0
end subroutine init_result_cat

subroutine independence_metropolis_sample(f,x0,proposal_mean,proposal_cov,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),proposal_mean(:),proposal_cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,info
   real(dp) :: x(size(x0)),prop(size(x0)),lp,lpp,lqx,lqp
   real(dp) :: inv(size(x0),size(x0)),d(size(x0)),ldet
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   call inverse_and_logdet(proposal_cov,inv,ldet,info)
   if(info/=0) return
   x=x0; lp=f(x); lqx=gaussian_logkernel(x,proposal_mean,inv,ldet); k=0
   do it=1,n_iter
      call rand_mvn(proposal_mean,proposal_cov,prop,info); if(info/=0) cycle
      lpp=f(prop); lqp=gaussian_logkernel(prop,proposal_mean,inv,ldet); res%proposed=res%proposed+1
      if(log(rand_uniform())<lpp-lp+lqx-lqp) then
         x=prop; lp=lpp; lqx=lqp; res%accepted=res%accepted+1
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine independence_metropolis_sample

subroutine multiple_try_metropolis_sample(f,x0,proposal_cov,ntry,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),proposal_cov(:,:)
   integer, intent(in) :: ntry,n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,m,j,sel,info
   real(dp) :: x(size(x0)),lp,step(size(x0)),sumw,sumback,u
   real(dp), allocatable :: cand(:,:),lpc(:),w(:),back(:,:),lpback(:)
   p=size(x0); m=max(2,ntry); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   allocate(cand(m,p),lpc(m),w(m),back(m,p),lpback(m)); x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,m
         call rand_mvn(0.0_dp*x,proposal_cov,step,info); cand(j,:)=x+step; lpc(j)=f(cand(j,:))
      end do
      call normalized_exp_weights(lpc,w); u=rand_uniform(); sel=weighted_index(w,u)
      back(1,:)=x; lpback(1)=lp
      do j=2,m
         call rand_mvn(0.0_dp*x,proposal_cov,step,info); back(j,:)=cand(sel,:)+step; lpback(j)=f(back(j,:))
      end do
      sumw=logsumexp(lpc); sumback=logsumexp(lpback); res%proposed=res%proposed+1
      if(log(rand_uniform())<min(0.0_dp,sumw-sumback)) then
         x=cand(sel,:); lp=lpc(sel); res%accepted=res%accepted+1
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine multiple_try_metropolis_sample

subroutine hit_and_run_sample(f,x0,step_sd,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),step_sd
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,j
   real(dp) :: x(size(x0)),prop(size(x0)),dir(size(x0)),normd,lp,lpp,step
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p; dir(j)=rand_normal(); end do
      normd=sqrt(max(dot_product(dir,dir),tiny(1.0_dp))); dir=dir/normd
      step=step_sd*rand_normal(); prop=x+step*dir; lpp=f(prop); res%proposed=res%proposed+1
      if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine hit_and_run_sample

subroutine charm_sample(f,x0,step_sd,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),step_sd(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,j,k,nkeep
   real(dp) :: x(size(x0)),prop(size(x0)),lp,lpp
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         prop=x; prop(j)=x(j)+step_sd(j)*rand_normal(); lpp=f(prop); res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      end do
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine charm_sample

subroutine adaptive_mwg_sample(f,x0,initial_sd,n_iter,burn,thin,res,target,adapt_every,final_sd)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),initial_sd(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   real(dp), intent(in), optional :: target
   integer, intent(in), optional :: adapt_every
   real(dp), allocatable, intent(out), optional :: final_sd(:)
   integer :: p,it,j,k,nkeep,period
   integer, allocatable :: aacc(:),aprop(:)
   real(dp) :: x(size(x0)),prop(size(x0)),lp,lpp,targ,eta
   real(dp), allocatable :: sd(:)
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   allocate(sd(p),aacc(p),aprop(p)); sd=max(initial_sd,1.0e-8_dp); aacc=0; aprop=0
   targ=0.44_dp; if(present(target)) targ=target
   period=25; if(present(adapt_every)) period=max(1,adapt_every)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         prop=x; prop(j)=x(j)+sd(j)*rand_normal(); lpp=f(prop); aprop(j)=aprop(j)+1; res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then
            x=prop; lp=lpp; aacc(j)=aacc(j)+1; res%accepted=res%accepted+1
         end if
      end do
      if(it<=burn .and. mod(it,period)==0) then
         eta=min(0.25_dp,1.0_dp/sqrt(real(max(1,it/period),dp)))
         do j=1,p
            if(aprop(j)>0) sd(j)=sd(j)*exp(eta*(real(aacc(j),dp)/real(aprop(j),dp)-targ))
         end do
         aacc=0; aprop=0
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_sd)) then; allocate(final_sd(p)); final_sd=sd; end if
end subroutine adaptive_mwg_sample

subroutine adaptive_mixture_metropolis_sample(f,x0,base_cov,n_iter,burn,thin,res,mix_prob,small_scale,adapt_start)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),base_cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   real(dp), intent(in), optional :: mix_prob,small_scale
   integer, intent(in), optional :: adapt_start
   integer :: p,it,k,nkeep,info,astart
   real(dp) :: x(size(x0)),prop(size(x0)),step(size(x0)),lp,lpp,pmix,ss,sc
   real(dp), allocatable :: cov(:,:),history(:,:)
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   pmix=0.05_dp; if(present(mix_prob)) pmix=min(max(mix_prob,0.0_dp),1.0_dp)
   ss=0.1_dp; if(present(small_scale)) ss=max(small_scale,1.0e-6_dp)
   astart=max(20,2*p); if(present(adapt_start)) astart=max(2,adapt_start)
   allocate(cov(p,p),history(n_iter,p)); cov=base_cov; call make_positive_definite(cov)
   x=x0; lp=f(x); k=0; sc=2.38_dp**2/real(p,dp)
   do it=1,n_iter
      if(it>astart .and. mod(it,20)==0) then
         call sample_covariance(history(1:it-1,:),cov); cov=sc*cov+1.0e-8_dp*identity_local(p)
         call make_positive_definite(cov)
      end if
      if(rand_uniform()<pmix) then
         call rand_mvn(0.0_dp*x,ss*ss*base_cov,step,info)
      else
         call rand_mvn(0.0_dp*x,cov,step,info)
      end if
      prop=x+step; lpp=f(prop); res%proposed=res%proposed+1
      if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      history(it,:)=x
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine adaptive_mixture_metropolis_sample

subroutine sgld_sample(f,x0,step0,n_iter,burn,thin,res,decay)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),step0
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   real(dp), intent(in), optional :: decay
   integer :: p,it,j,k,nkeep
   real(dp) :: x(size(x0)),g(size(x0)),eps,alpha,lp
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   alpha=0.55_dp; if(present(decay)) alpha=decay
   x=x0; k=0
   do it=1,n_iter
      eps=step0/(real(it,dp)**alpha); call numerical_gradient(f,x,g)
      do j=1,p; x(j)=x(j)+0.5_dp*eps*g(j)+sqrt(eps)*rand_normal(); end do
      lp=f(x); res%proposed=res%proposed+1; res%accepted=res%accepted+1
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   res%acceptance_rate=1.0_dp
end subroutine sgld_sample

subroutine griddy_gibbs_sample(f,x0,width,grid_points,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),width(:)
   integer, intent(in) :: grid_points,n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,j,k,nkeep,m,g,sel
   real(dp) :: x(size(x0)),trial(size(x0)),lp,u
   real(dp), allocatable :: vals(:),logw(:),w(:)
   p=size(x0); m=max(5,grid_points); if(mod(m,2)==0) m=m+1
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   allocate(vals(m),logw(m),w(m)); x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         do g=1,m
            vals(g)=x(j)+width(j)*(2.0_dp*real(g-1,dp)/real(m-1,dp)-1.0_dp)
            trial=x; trial(j)=vals(g); logw(g)=f(trial)
         end do
         call normalized_exp_weights(logw,w); u=rand_uniform(); sel=weighted_index(w,u)
         x(j)=vals(sel); lp=logw(sel); res%proposed=res%proposed+m; res%accepted=res%accepted+1
      end do
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine griddy_gibbs_sample

subroutine adaptive_griddy_gibbs_sample(f,x0,width,grid_points,n_iter,burn,thin,res,adapt_every,final_width)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),width(:)
   integer, intent(in) :: grid_points,n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: adapt_every
   real(dp), allocatable, intent(out), optional :: final_width(:)
   integer :: p,it,j,k,nkeep,m,g,sel,period
   real(dp) :: x(size(x0)),trial(size(x0)),lp,u,meanx
   real(dp), allocatable :: vals(:),logw(:),w(:),wd(:),hist(:,:)
   p=size(x0); m=max(5,grid_points); if(mod(m,2)==0) m=m+1
   period=50; if(present(adapt_every)) period=max(5,adapt_every)
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   allocate(vals(m),logw(m),w(m),wd(p),hist(max(1,burn),p)); wd=max(width,1.0e-6_dp)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         do g=1,m
            vals(g)=x(j)+wd(j)*(2.0_dp*real(g-1,dp)/real(m-1,dp)-1.0_dp)
            trial=x; trial(j)=vals(g); logw(g)=f(trial)
         end do
         call normalized_exp_weights(logw,w); u=rand_uniform(); sel=weighted_index(w,u)
         x(j)=vals(sel); lp=logw(sel); res%proposed=res%proposed+m; res%accepted=res%accepted+1
      end do
      if(it<=burn) hist(it,:)=x
      if(it<=burn .and. mod(it,period)==0 .and. it>5) then
         do j=1,p
            meanx=sum(hist(1:it,j))/real(it,dp)
            wd(j)=max(width(j)*0.1_dp,4.0_dp*sqrt(sum((hist(1:it,j)-meanx)**2)/real(it-1,dp)))
         end do
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_width)) then; allocate(final_width(p)); final_width=wd; end if
end subroutine adaptive_griddy_gibbs_sample

subroutine ohss_sample(f,x0,width,n_iter,burn,thin,res,max_steps)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),width
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: max_steps
   integer :: p,it,k,nkeep,j,m,steps
   real(dp) :: x(size(x0)),trial(size(x0)),dir(size(x0)),lp,logy,l,r,t,normd
   p=size(x0); m=100; if(present(max_steps)) m=max(1,max_steps)
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p; dir(j)=rand_normal(); end do; normd=sqrt(dot_product(dir,dir)); dir=dir/max(normd,tiny(1.0_dp))
      logy=lp+log(rand_uniform()); l=-width*rand_uniform(); r=l+width
      steps=0
      do while(steps<m); trial=x+l*dir; if(f(trial)<=logy) exit; l=l-width; steps=steps+1; end do
      steps=0
      do while(steps<m); trial=x+r*dir; if(f(trial)<=logy) exit; r=r+width; steps=steps+1; end do
      do
         t=l+(r-l)*rand_uniform(); trial=x+t*dir; res%proposed=res%proposed+1
         if(f(trial)>=logy) then; x=trial; lp=f(x); res%accepted=res%accepted+1; exit; end if
         if(t<0.0_dp) then; l=t; else; r=t; end if
      end do
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine ohss_sample

subroutine uess_sample(f,x0,basis_cov,width,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),basis_cov(:,:),width
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,j,info,col
   real(dp) :: x(size(x0)),trial(size(x0)),dir(size(x0)),lp,logy,l,r,t,normd
   real(dp), allocatable :: lmat(:,:)
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   allocate(lmat(p,p)); call chol_lower(basis_cov,lmat,info); if(info/=0) lmat=identity_local(p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      col=1+mod(it-1,p); dir=lmat(:,col); normd=sqrt(dot_product(dir,dir)); dir=dir/max(normd,tiny(1.0_dp))
      logy=lp+log(rand_uniform()); l=-width*rand_uniform(); r=l+width
      do j=1,50; trial=x+l*dir; if(f(trial)<=logy) exit; l=l-width; end do
      do j=1,50; trial=x+r*dir; if(f(trial)<=logy) exit; r=r+width; end do
      do
         t=l+(r-l)*rand_uniform(); trial=x+t*dir; res%proposed=res%proposed+1
         if(f(trial)>=logy) then; x=trial; lp=f(x); res%accepted=res%accepted+1; exit; end if
         if(t<0.0_dp) then; l=t; else; r=t; end if
      end do
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine uess_sample

subroutine afss_sample(f,x0,initial_cov,width,n_iter,burn,thin,res,adapt_every)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),initial_cov(:,:),width
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: adapt_every
   integer :: p,it,k,nkeep,period,info,col,j
   real(dp) :: x(size(x0)),trial(size(x0)),dir(size(x0)),lp,logy,l,r,t,normd
   real(dp), allocatable :: cov(:,:),fac(:,:),hist(:,:)
   p=size(x0); period=50; if(present(adapt_every)) period=max(10,adapt_every)
   nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   allocate(cov(p,p),fac(p,p),hist(max(1,burn),p)); cov=initial_cov; call make_positive_definite(cov)
   call chol_lower(cov,fac,info); x=x0; lp=f(x); k=0
   do it=1,n_iter
      col=1+mod(it-1,p); dir=fac(:,col); normd=sqrt(dot_product(dir,dir)); dir=dir/max(normd,tiny(1.0_dp))
      logy=lp+log(rand_uniform()); l=-width*rand_uniform(); r=l+width
      do j=1,50; trial=x+l*dir; if(f(trial)<=logy) exit; l=l-width; end do
      do j=1,50; trial=x+r*dir; if(f(trial)<=logy) exit; r=r+width; end do
      do
         t=l+(r-l)*rand_uniform(); trial=x+t*dir; res%proposed=res%proposed+1
         if(f(trial)>=logy) then; x=trial; lp=f(x); res%accepted=res%accepted+1; exit; end if
         if(t<0.0_dp) then; l=t; else; r=t; end if
      end do
      if(it<=burn) hist(it,:)=x
      if(it<=burn .and. it>p+2 .and. mod(it,period)==0) then
         call sample_covariance(hist(1:it,:),cov); cov=cov+1.0e-8_dp*identity_local(p); call make_positive_definite(cov)
         call chol_lower(cov,fac,info)
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine afss_sample

subroutine mcmcmc_sample(f,chains,temperatures,n_iter,burn,thin,out,swap_rate,proposal_sd,swap_every)
   procedure(log_target_iface) :: f
   real(dp), intent(inout) :: chains(:,:)
   real(dp), intent(in) :: temperatures(:),proposal_sd
   integer, intent(in) :: n_iter,burn,thin
   real(dp), allocatable, intent(out) :: out(:,:,:)
   real(dp), intent(out) :: swap_rate
   integer, intent(in), optional :: swap_every
   integer :: nc,p,it,c,j,k,nkeep,se,swaps,swapacc
   real(dp) :: lp(size(chains,1)),prop(size(chains,2)),lpp,loga,tmpv(size(chains,2)),u
   nc=size(chains,1); p=size(chains,2); se=10; if(present(swap_every)) se=max(1,swap_every)
   nkeep=max(0,(n_iter-burn)/max(1,thin)); allocate(out(nkeep,nc,p)); k=0; swaps=0; swapacc=0
   do c=1,nc; lp(c)=f(chains(c,:)); end do
   do it=1,n_iter
      do c=1,nc
         prop=chains(c,:); do j=1,p; prop(j)=prop(j)+proposal_sd*rand_normal(); end do
         lpp=f(prop); loga=(lpp-lp(c))/temperatures(c)
         if(log(rand_uniform())<loga) then; chains(c,:)=prop; lp(c)=lpp; end if
      end do
      if(mod(it,se)==0 .and. nc>1) then
         c=1+int(rand_uniform()*real(nc-1,dp)); swaps=swaps+1
         loga=(1.0_dp/temperatures(c)-1.0_dp/temperatures(c+1))*(lp(c+1)-lp(c))
         if(log(rand_uniform())<loga) then
            tmpv=chains(c,:); chains(c,:)=chains(c+1,:); chains(c+1,:)=tmpv
            u=lp(c); lp(c)=lp(c+1); lp(c+1)=u; swapacc=swapacc+1
         end if
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; out(k,:,:)=chains
      end if
   end do
   swap_rate=real(swapacc,dp)/real(max(1,swaps),dp)
end subroutine mcmcmc_sample

subroutine random_dive_sample(f,x0,scale,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),scale
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,j,k,nkeep,m,tmp,ord(size(x0)),sgn(size(x0))
   real(dp) :: x(size(x0)),prop(size(x0)),lp,lpp,u(size(x0)),eps(size(x0)),loga,jac,r
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   x=x0; lp=f(x); k=0
   ! scale is retained for source compatibility with the initial v0.3 draft.
   ! Upstream LaplacesDemon RDMH has no tuning scale.
   r=scale*0.0_dp
   do it=1,n_iter
      do j=1,p
         sgn(j)=merge(1,-1,rand_uniform()<0.5_dp)
         u(j)=2.0_dp*rand_uniform()-1.0_dp
         if(abs(u(j))<epsilon(1.0_dp)) u(j)=sign(epsilon(1.0_dp),u(j)+epsilon(1.0_dp))
         if(sgn(j)==1) then
            eps(j)=u(j)
         else
            eps(j)=1.0_dp/u(j)
         end if
         ord(j)=j
      end do
      do j=p,2,-1
         m=1+int(rand_uniform()*real(j,dp))
         tmp=ord(j); ord(j)=ord(m); ord(m)=tmp
      end do
      do m=1,p
         j=ord(m); prop=x; prop(j)=x(j)*eps(j); lpp=f(prop)
         if(x(j)/=0.0_dp .and. prop(j)/=0.0_dp) then
            jac=log(abs(prop(j)/x(j)))
            if(jac/=jac .or. abs(jac)>=huge(jac)) jac=0.0_dp
         else
            jac=0.0_dp
         end if
         loga=lpp-lp+jac; res%proposed=res%proposed+1
         if(log(rand_uniform())<loga) then
            x=prop; lp=lpp; res%accepted=res%accepted+1
         end if
      end do
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine random_dive_sample

subroutine rss_sample(f,x0,width,lower,upper,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),width(:),lower(:),upper(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,j,k,nkeep
   real(dp) :: x(size(x0)),trial(size(x0)),lp,logy,l,r,t
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         logy=lp+log(rand_uniform()); l=max(lower(j),x(j)-width(j)*rand_uniform()); r=min(upper(j),l+width(j))
         do
            t=l+(r-l)*rand_uniform(); trial=x; trial(j)=reflect_box(t,lower(j),upper(j)); res%proposed=res%proposed+1
            if(f(trial)>=logy) then; x=trial; lp=f(x); res%accepted=res%accepted+1; exit; end if
            if(trial(j)<x(j)) then; l=t; else; r=t; end if
            if(abs(r-l)<1.0e-12_dp) exit
         end do
      end do
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine rss_sample

subroutine tempered_hmc_sample(f,x0,step_size,n_leapfrog,temperature,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),step_size,temperature
   integer, intent(in) :: n_leapfrog,n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,l,j
   real(dp) :: x(size(x0)),q(size(x0)),mom(size(x0)),mom0(size(x0)),g(size(x0)),lp,lpp,h0,h1,temp
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_cat(res,nkeep,p)
   temp=max(temperature,1.0_dp); x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p; mom0(j)=sqrt(temp)*rand_normal(); end do
      mom=mom0; q=x; call numerical_gradient(f,q,g); mom=mom+0.5_dp*step_size*g
      do l=1,n_leapfrog
         q=q+step_size*mom/temp; call numerical_gradient(f,q,g)
         if(l<n_leapfrog) mom=mom+step_size*g
      end do
      mom=mom+0.5_dp*step_size*g; lpp=f(q)
      h0=-lp+0.5_dp*dot_product(mom0,mom0)/temp; h1=-lpp+0.5_dp*dot_product(mom,mom)/temp
      res%proposed=res%proposed+1
      if(log(rand_uniform())<min(0.0_dp,h0-h1)) then; x=q; lp=lpp; res%accepted=res%accepted+1; end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine tempered_hmc_sample

subroutine inverse_and_logdet(a,ainv,ldet,info)
   use ld_linalg, only: inverse_spd, logdet_spd
   real(dp), intent(in) :: a(:,:)
   real(dp), intent(out) :: ainv(:,:),ldet
   integer, intent(out) :: info
   call inverse_spd(a,ainv,info); if(info==0) ldet=logdet_spd(a,info)
end subroutine inverse_and_logdet

pure function gaussian_logkernel(x,mu,inv,ldet) result(v)
   real(dp), intent(in) :: x(:),mu(:),inv(:,:),ldet
   real(dp) :: v,d(size(x))
   d=x-mu; v=-0.5_dp*(real(size(x),dp)*log(2.0_dp*pi)+ldet+dot_product(d,matmul(inv,d)))
end function gaussian_logkernel

subroutine normalized_exp_weights(logw,w)
   real(dp), intent(in) :: logw(:)
   real(dp), intent(out) :: w(:)
   real(dp) :: m,s
   m=maxval(logw); w=exp(logw-m); s=sum(w)
   if(s<=0.0_dp .or. s/=s) then; w=1.0_dp/real(size(w),dp); else; w=w/s; end if
end subroutine normalized_exp_weights

integer function weighted_index(w,u) result(idx)
   real(dp), intent(in) :: w(:),u
   real(dp) :: c
   integer :: i
   c=0.0_dp; idx=size(w)
   do i=1,size(w); c=c+w(i); if(u<=c) then; idx=i; return; end if; end do
end function weighted_index

pure function logsumexp(x) result(v)
   real(dp), intent(in) :: x(:)
   real(dp) :: v,m
   m=maxval(x); v=m+log(sum(exp(x-m)))
end function logsumexp

pure function identity_local(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
end function identity_local

pure function reflect_box(x,lo,hi) result(y)
   real(dp), intent(in) :: x,lo,hi
   real(dp) :: y,w,z
   if(hi<=lo) then; y=lo; return; end if
   w=hi-lo; z=modulo(x-lo,2.0_dp*w)
   if(z<=w) then; y=lo+z; else; y=hi-(z-w); end if
end function reflect_box

end module ld_mcmc_catalog
