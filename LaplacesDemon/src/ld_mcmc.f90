module ld_mcmc
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface
use ld_random, only: rand_uniform, rand_normal, rand_mvn
use ld_linalg, only: sample_covariance, make_positive_definite
use ld_numerics, only: numerical_gradient
implicit none
private
public :: mcmc_result_t, rwm_sample, mwg_sample, adaptive_metropolis_sample
public :: mala_sample, hmc_sample, slice_sample, elliptical_slice_sample
public :: aies_sample, demc_sample

type :: mcmc_result_t
   real(dp), allocatable :: chain(:,:), logp(:)
   real(dp) :: acceptance_rate=0.0_dp
   integer :: proposed=0, accepted=0
end type
contains
subroutine init_result(res,nkeep,p)
   type(mcmc_result_t),intent(out)::res
   integer,intent(in)::nkeep,p
   allocate(res%chain(nkeep,p),res%logp(nkeep)); res%chain=0.0_dp; res%logp=0.0_dp
end subroutine init_result

subroutine rwm_sample(f,x0,proposal_cov,n_iter,burn,thin,res)
   procedure(log_target_iface)::f
   real(dp),intent(in)::x0(:),proposal_cov(:,:)
   integer,intent(in)::n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer::p,it,k,nkeep,info
   real(dp)::x(size(x0)),prop(size(x0)),step(size(x0)),lp,lpp
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      call rand_mvn(0.0_dp*x,proposal_cov,step,info); prop=x+step; lpp=f(prop); res%proposed=res%proposed+1
      if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine rwm_sample

subroutine mwg_sample(f,x0,proposal_sd,n_iter,burn,thin,res)
   procedure(log_target_iface)::f
   real(dp),intent(in)::x0(:),proposal_sd(:)
   integer,intent(in)::n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer::p,it,j,k,nkeep
   real(dp)::x(size(x0)),prop(size(x0)),lp,lpp
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p); x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         prop=x; prop(j)=x(j)+proposal_sd(j)*rand_normal(); lpp=f(prop); res%proposed=res%proposed+1
         if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      end do
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine mwg_sample

subroutine adaptive_metropolis_sample(f,x0,initial_cov,n_iter,burn,thin,res,adapt_start,scale)
   procedure(log_target_iface)::f
   real(dp),intent(in)::x0(:),initial_cov(:,:)
   integer,intent(in)::n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer,intent(in),optional::adapt_start
   real(dp),intent(in),optional::scale
   integer::p,it,k,nkeep,astart,info
   real(dp)::x(size(x0)),prop(size(x0)),step(size(x0)),lp,lpp,sc
   real(dp),allocatable::history(:,:),cov(:,:)
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p)
   astart=max(20,2*p); if(present(adapt_start)) astart=adapt_start
   sc=2.38_dp**2/real(p,dp); if(present(scale)) sc=scale
   allocate(history(n_iter,p),cov(p,p)); cov=initial_cov; call make_positive_definite(cov)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      if(it>astart .and. mod(it,10)==0) then
         call sample_covariance(history(1:it-1,:),cov); cov=sc*cov
         cov=cov+1e-10_dp*identity(p); call make_positive_definite(cov)
      end if
      call rand_mvn(0.0_dp*x,cov,step,info); prop=x+step; lpp=f(prop); res%proposed=res%proposed+1
      if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      history(it,:)=x
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
contains
   pure function identity(n) result(a)
      integer,intent(in)::n; real(dp)::a(n,n); integer::i
      a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
   end function identity
end subroutine adaptive_metropolis_sample

subroutine mala_sample(f,x0,step_size,n_iter,burn,thin,res)
   procedure(log_target_iface)::f
   real(dp),intent(in)::x0(:),step_size
   integer,intent(in)::n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer::p,it,k,nkeep,j
   real(dp)::x(size(x0)),prop(size(x0)),g(size(x0)),gp(size(x0)),meanf(size(x0)),meanr(size(x0)),lp,lpp,loga,s2
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p); s2=step_size*step_size
   x=x0; lp=f(x); call numerical_gradient(f,x,g); k=0
   do it=1,n_iter
      meanf=x+0.5_dp*s2*g; do j=1,p; prop(j)=meanf(j)+step_size*rand_normal(); end do
      lpp=f(prop); call numerical_gradient(f,prop,gp); meanr=prop+0.5_dp*s2*gp
      loga=lpp-lp-0.5_dp*sum((x-meanr)**2)/s2+0.5_dp*sum((prop-meanf)**2)/s2; res%proposed=res%proposed+1
      if(log(rand_uniform())<loga) then; x=prop; lp=lpp; g=gp; res%accepted=res%accepted+1; end if
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine mala_sample

subroutine hmc_sample(f,x0,step_size,n_leapfrog,n_iter,burn,thin,res)
   procedure(log_target_iface)::f
   real(dp),intent(in)::x0(:),step_size
   integer,intent(in)::n_leapfrog,n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer::p,it,k,nkeep,l,j
   real(dp)::x(size(x0)),q(size(x0)),mom(size(x0)),mom0(size(x0)),g(size(x0)),lp,lpp,h0,h1
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p); x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p; mom0(j)=rand_normal(); end do; mom=mom0; q=x; call numerical_gradient(f,q,g); mom=mom+0.5_dp*step_size*g
      do l=1,n_leapfrog
         q=q+step_size*mom; call numerical_gradient(f,q,g); if(l<n_leapfrog) mom=mom+step_size*g
      end do
      mom=mom+0.5_dp*step_size*g; lpp=f(q); h0=-lp+0.5_dp*sum(mom0*mom0); h1=-lpp+0.5_dp*sum(mom*mom); res%proposed=res%proposed+1
      if(log(rand_uniform())<min(0.0_dp,h0-h1)) then; x=q; lp=lpp; res%accepted=res%accepted+1; end if
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine hmc_sample

subroutine slice_sample(f,x0,width,n_iter,burn,thin,res,max_steps)
   procedure(log_target_iface)::f
   real(dp),intent(in)::x0(:),width(:)
   integer,intent(in)::n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer,intent(in),optional::max_steps
   integer::p,it,j,k,nkeep,m,steps
   real(dp)::x(size(x0)),trial(size(x0)),lp,logy,l,r,u
   p=size(x0); m=100; if(present(max_steps)) m=max_steps; nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do j=1,p
         logy=lp+log(rand_uniform()); u=rand_uniform(); l=x(j)-u*width(j); r=l+width(j)
         trial=x; steps=0
         do while(steps<m); trial(j)=l; if(f(trial)<=logy) exit; l=l-width(j); steps=steps+1; end do
         steps=0; do while(steps<m); trial(j)=r; if(f(trial)<=logy) exit; r=r+width(j); steps=steps+1; end do
         do
            trial=x; trial(j)=l+(r-l)*rand_uniform(); res%proposed=res%proposed+1
            if(f(trial)>=logy) then; x=trial; lp=f(x); res%accepted=res%accepted+1; exit; end if
            if(trial(j)<x(j)) then; l=trial(j); else; r=trial(j); end if
         end do
      end do
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=lp; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine slice_sample

subroutine elliptical_slice_sample(loglik,x0,prior_cov,n_iter,burn,thin,res)
   procedure(log_target_iface)::loglik
   real(dp),intent(in)::x0(:),prior_cov(:,:)
   integer,intent(in)::n_iter,burn,thin
   type(mcmc_result_t),intent(out)::res
   integer::p,it,k,nkeep,info
   real(dp)::x(size(x0)),nu(size(x0)),prop(size(x0)),ll,logy,theta,tmin,tmax
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(thin,1)); call init_result(res,nkeep,p); x=x0; ll=loglik(x); k=0
   do it=1,n_iter
      call rand_mvn(0.0_dp*x,prior_cov,nu,info); logy=ll+log(rand_uniform()); theta=2.0_dp*acos(-1.0_dp)*rand_uniform()
      tmin=theta-2.0_dp*acos(-1.0_dp); tmax=theta
      do
         prop=x*cos(theta)+nu*sin(theta); res%proposed=res%proposed+1
         if(loglik(prop)>logy) then; x=prop; ll=loglik(x); res%accepted=res%accepted+1; exit; end if
         if(theta<0.0_dp) then; tmin=theta; else; tmax=theta; end if; theta=tmin+(tmax-tmin)*rand_uniform()
      end do
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; res%chain(k,:)=x; res%logp(k)=ll; end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine elliptical_slice_sample

subroutine aies_sample(f,walkers,n_iter,burn,thin,chain,acceptance,a_scale)
   procedure(log_target_iface)::f
   real(dp),intent(inout)::walkers(:,:)
   integer,intent(in)::n_iter,burn,thin
   real(dp),allocatable,intent(out)::chain(:,:,:)
   real(dp),intent(out)::acceptance
   real(dp),intent(in),optional::a_scale
   integer::nw,p,it,i,j,k,nkeep,acc,propn
   real(dp)::a,z,lp(size(walkers,1)),lpp,proposal(size(walkers,2)),u
   nw=size(walkers,1); p=size(walkers,2); a=2.0_dp; if(present(a_scale)) a=a_scale
   nkeep=max(0,(n_iter-burn)/max(thin,1)); allocate(chain(nkeep,nw,p)); k=0; acc=0; propn=0
   do i=1,nw; lp(i)=f(walkers(i,:)); end do
   do it=1,n_iter
      do i=1,nw
         j=1+int(rand_uniform()*real(nw-1,dp)); if(j>=i) j=j+1
         u=rand_uniform(); z=((a-1.0_dp)*u+1.0_dp)**2/a
         proposal=walkers(j,:)+z*(walkers(i,:)-walkers(j,:)); lpp=f(proposal); propn=propn+1
         if(log(rand_uniform())<(p-1)*log(z)+lpp-lp(i)) then; walkers(i,:)=proposal; lp(i)=lpp; acc=acc+1; end if
      end do
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; chain(k,:,:)=walkers; end if
   end do
   acceptance=real(acc,dp)/real(max(1,propn),dp)
end subroutine aies_sample

subroutine demc_sample(f,chains,n_iter,burn,thin,out,acceptance,gamma_scale,noise_sd)
   procedure(log_target_iface)::f
   real(dp),intent(inout)::chains(:,:)
   integer,intent(in)::n_iter,burn,thin
   real(dp),allocatable,intent(out)::out(:,:,:)
   real(dp),intent(out)::acceptance
   real(dp),intent(in),optional::gamma_scale,noise_sd
   integer::nc,p,it,i,r1,r2,j,k,nkeep,acc,propn
   real(dp)::gamma,eps,lp(size(chains,1)),lpp,proposal(size(chains,2))
   nc=size(chains,1)
   p=size(chains,2)
   gamma=2.38_dp/sqrt(2.0_dp*real(p,dp))
   if(present(gamma_scale)) gamma=gamma_scale
   eps=1e-6_dp
   if(present(noise_sd)) eps=noise_sd
   nkeep=max(0,(n_iter-burn)/max(thin,1))
   allocate(out(nkeep,nc,p))
   k=0; acc=0; propn=0
   do i=1,nc; lp(i)=f(chains(i,:)); end do
   do it=1,n_iter
      do i=1,nc
         do; r1=1+int(rand_uniform()*real(nc,dp)); if(r1/=i) exit; end do
         do; r2=1+int(rand_uniform()*real(nc,dp)); if(r2/=i .and. r2/=r1) exit; end do
         proposal=chains(i,:)+gamma*(chains(r1,:)-chains(r2,:)); do j=1,p; proposal(j)=proposal(j)+eps*rand_normal(); end do
         lpp=f(proposal); propn=propn+1; if(log(rand_uniform())<lpp-lp(i)) then; chains(i,:)=proposal; lp(i)=lpp; acc=acc+1; end if
      end do
      if(it>burn .and. mod(it-burn,thin)==0 .and. k<nkeep) then; k=k+1; out(k,:,:)=chains; end if
   end do
   acceptance=real(acc,dp)/real(max(1,propn),dp)
end subroutine demc_sample
end module ld_mcmc
