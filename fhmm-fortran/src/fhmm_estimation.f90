! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_estimation
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fhmm_kinds, only: dp
   use fhmm_types, only: hmm_parameters, hhmm_parameters, fit_options
   use fhmm_types, only: hmm_fit_result, hhmm_fit_result
   use fhmm_types, only: dist_gamma, dist_poisson, dist_student_t
   use fhmm_parameters, only: pack_hmm_parameters,unpack_hmm_parameters,hmm_parameter_count
   use fhmm_parameters, only: pack_hhmm_parameters,unpack_hhmm_parameters,hhmm_parameter_count
   use fhmm_algorithms, only: hmm_log_likelihood,forward_backward,viterbi_decode
   use fhmm_hierarchical, only: hhmm_log_likelihood
   use fhmm_math, only: seed_rng,numerical_gradient,numerical_hessian,invert_matrix
   use fhmm_optimize, only: nelder_mead
   implicit none
   private
   public :: initial_hmm_parameters, initial_hhmm_parameters, fit_hmm, fit_hhmm

   real(dp),allocatable,save::ctx_obs(:),ctx_coarse(:),ctx_fine(:,:)
   integer,allocatable,save::ctx_lengths(:)
   integer,save::ctx_states=0,ctx_family=0,ctx_m=0,ctx_n=0,ctx_coarse_family=0,ctx_fine_family=0

contains

   function initial_hmm_parameters(observations,nstates,family) result(par)
      real(dp),intent(in)::observations(:)
      integer,intent(in)::nstates,family
      type(hmm_parameters)::par
      real(dp),allocatable::sorted(:)
      real(dp)::overall_sd
      integer::i,k,Tn
      Tn=size(observations);allocate(sorted(Tn));sorted=observations;call sort_real(sorted)
      par%distribution=family;allocate(par%gamma(nstates,nstates),par%mu(nstates),par%sigma(nstates),par%df(nstates))
      if(nstates==1)then
         par%gamma=1.0_dp
      else
         par%gamma=0.1_dp/real(nstates-1,dp)
         do i=1,nstates;par%gamma(i,i)=0.9_dp;end do
      end if
      overall_sd=sample_sd(observations);if(overall_sd<=1.0e-8_dp)overall_sd=1.0_dp
      do i=1,nstates
         k=max(1,min(Tn,nint((real(i,dp)-0.5_dp)*real(Tn,dp)/real(nstates,dp))))
         par%mu(i)=sorted(k)
      end do
      if(family==dist_gamma.or.family==dist_poisson)par%mu=max(par%mu,0.1_dp)
      par%sigma=max(overall_sd/real(max(1,nstates),dp),0.1_dp*overall_sd)
      par%df=8.0_dp
   end function initial_hmm_parameters

   function initial_hhmm_parameters(coarse_obs,fine_obs,lengths,m,n,coarse_family,fine_family) result(par)
      real(dp),intent(in)::coarse_obs(:),fine_obs(:,:)
      integer,intent(in)::lengths(:),m,n,coarse_family,fine_family
      type(hhmm_parameters)::par
      real(dp),allocatable::pooled(:)
      integer::t,k,total,s
      par%coarse=initial_hmm_parameters(coarse_obs,m,coarse_family);allocate(par%fine(m))
      total=sum(lengths);allocate(pooled(total));k=0
      do t=1,size(lengths);pooled(k+1:k+lengths(t))=fine_obs(t,1:lengths(t));k=k+lengths(t);end do
      do s=1,m;par%fine(s)=initial_hmm_parameters(pooled,n,fine_family);end do
   end function initial_hhmm_parameters

   function fit_hmm(observations,nstates,family,options,initial) result(out)
      real(dp),intent(in)::observations(:)
      integer,intent(in)::nstates,family
      type(fit_options),intent(in),optional::options
      type(hmm_parameters),intent(in),optional::initial
      type(hmm_fit_result)::out
      type(fit_options)::opt
      type(hmm_parameters)::start
      real(dp),allocatable::u0(:),u(:),best_u(:),g(:),hess(:,:),cov(:,:),trial(:)
      real(dp)::fval,best_f
      integer::run,it,best_it,p
      logical::ok,inv_ok
      opt=fit_options();if(present(options))opt=options
      if(size(observations)<2.or.nstates<1)then;out%message='invalid observations or state count';return;end if
      if(present(initial))then;start=initial;else;start=initial_hmm_parameters(observations,nstates,family);end if
      u0=pack_hmm_parameters(start);p=size(u0);allocate(best_u(p),u(p),trial(p));best_f=huge(1.0_dp)
      ctx_obs=observations;ctx_states=nstates;ctx_family=family;best_it=0;call seed_rng(opt%seed)
      do run=1,max(1,opt%runs)
         trial=u0
         if(run>1)call jitter_vector(trial,opt%initial_jitter)
         call nelder_mead(hmm_objective,trial,opt%max_iterations,opt%x_tolerance,opt%f_tolerance,u,fval,it,ok)
         if(fval<best_f)then;best_f=fval;best_u=u;best_it=it;end if
      end do
      if(.not.ieee_is_finite(best_f).or.best_f>=huge(1.0_dp)/1000.0_dp)then;out%message='optimization failed';return;end if
      out%parameters=unpack_hmm_parameters(best_u,nstates,family,.true.);out%unconstrained=best_u
      out%log_likelihood=-best_f;out%iterations=best_it;out%inference=forward_backward(observations,out%parameters)
      out%decoding=viterbi_decode(observations,out%parameters)
      allocate(g(p))
      call numerical_gradient(hmm_objective,best_u,g)
      out%gradient=g
      if(opt%compute_hessian)then
         allocate(hess(p,p),cov(p,p),out%standard_error(p));call numerical_hessian(hmm_objective,best_u,hess)
         call invert_matrix(hess,cov,inv_ok)
         if(inv_ok)then
            out%covariance=cov;out%standard_error=sqrt(max(diagonal(cov),0.0_dp))
         else
            allocate(out%covariance(p,p));out%covariance=0.0_dp;out%standard_error=0.0_dp
         end if
      end if
      out%aic=-2.0_dp*out%log_likelihood+2.0_dp*real(p,dp)
      out%bic=-2.0_dp*out%log_likelihood+real(p,dp)*log(real(size(observations),dp))
      out%ok=.true.;out%message='ok'
   end function fit_hmm

   function fit_hhmm(coarse_obs,fine_obs,lengths,m,n,coarse_family,fine_family,options,initial) result(out)
      real(dp),intent(in)::coarse_obs(:),fine_obs(:,:)
      integer,intent(in)::lengths(:),m,n,coarse_family,fine_family
      type(fit_options),intent(in),optional::options
      type(hhmm_parameters),intent(in),optional::initial
      type(hhmm_fit_result)::out
      type(fit_options)::opt
      type(hhmm_parameters)::start
      real(dp),allocatable::u0(:),u(:),best_u(:),trial(:),hess(:,:),cov(:,:)
      real(dp)::fval,best_f
      integer::run,it,best_it,p
      logical::ok,inv_ok
      opt=fit_options();if(present(options))opt=options
      if(present(initial))then
         start=initial
      else
         start=initial_hhmm_parameters(coarse_obs,fine_obs,lengths,m,n,coarse_family,fine_family)
      end if
      u0=pack_hhmm_parameters(start);p=size(u0);allocate(u(p),best_u(p),trial(p));best_f=huge(1.0_dp)
      ctx_coarse=coarse_obs;ctx_fine=fine_obs;ctx_lengths=lengths;ctx_m=m;ctx_n=n
      ctx_coarse_family=coarse_family;ctx_fine_family=fine_family;best_it=0;call seed_rng(opt%seed)
      do run=1,max(1,opt%runs)
         trial=u0;if(run>1)call jitter_vector(trial,opt%initial_jitter)
         call nelder_mead(hhmm_objective,trial,opt%max_iterations,opt%x_tolerance,opt%f_tolerance,u,fval,it,ok)
         if(fval<best_f)then;best_f=fval;best_u=u;best_it=it;end if
      end do
      if(.not.ieee_is_finite(best_f))then;out%message='optimization failed';return;end if
      out%parameters=unpack_hhmm_parameters(best_u,m,n,coarse_family,fine_family,.true.);out%unconstrained=best_u
      out%log_likelihood=-best_f;out%iterations=best_it
      if(opt%compute_hessian)then
         allocate(hess(p,p),cov(p,p),out%standard_error(p));call numerical_hessian(hhmm_objective,best_u,hess)
         call invert_matrix(hess,cov,inv_ok)
         if(inv_ok)then;out%covariance=cov;out%standard_error=sqrt(max(diagonal(cov),0.0_dp))
         else;allocate(out%covariance(p,p));out%covariance=0.0_dp;out%standard_error=0.0_dp;end if
      end if
      out%aic=-2.0_dp*out%log_likelihood+2.0_dp*real(p,dp)
      out%bic=-2.0_dp*out%log_likelihood+real(p,dp)*log(real(size(coarse_obs)+sum(lengths),dp))
      out%ok=.true.;out%message='ok'
   end function fit_hhmm

   real(dp) function hmm_objective(u) result(value)
      real(dp),intent(in)::u(:)
      type(hmm_parameters)::par
      par=unpack_hmm_parameters(u,ctx_states,ctx_family,.true.)
      value=-hmm_log_likelihood(ctx_obs,par)
      if(.not.ieee_is_finite(value))value=huge(1.0_dp)/100.0_dp
   end function hmm_objective

   real(dp) function hhmm_objective(u) result(value)
      real(dp),intent(in)::u(:)
      type(hhmm_parameters)::par
      par=unpack_hhmm_parameters(u,ctx_m,ctx_n,ctx_coarse_family,ctx_fine_family,.true.)
      value=-hhmm_log_likelihood(ctx_coarse,ctx_fine,ctx_lengths,par)
      if(.not.ieee_is_finite(value))value=huge(1.0_dp)/100.0_dp
   end function hhmm_objective

   subroutine jitter_vector(x,factor)
      real(dp),intent(inout)::x(:)
      real(dp),intent(in)::factor
      real(dp)::u
      integer::i
      do i=1,size(x);call random_number(u);x(i)=x(i)+factor*(2.0_dp*u-1.0_dp)*max(1.0_dp,abs(x(i)));end do
   end subroutine jitter_vector

   subroutine sort_real(x)
      real(dp),intent(inout)::x(:)
      real(dp)::v
      integer::i,j
      do i=2,size(x);v=x(i);j=i-1;do while(j>=1);if(x(j)<=v)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=v;end do
   end subroutine sort_real

   pure real(dp) function sample_sd(x) result(s)
      real(dp),intent(in)::x(:)
      real(dp)::m
      if(size(x)<2)then;s=0.0_dp;else;m=sum(x)/real(size(x),dp);s=sqrt(sum((x-m)**2)/real(size(x)-1,dp));end if
   end function sample_sd

   pure function diagonal(a) result(d)
      real(dp),intent(in)::a(:,:)
      real(dp)::d(min(size(a,1),size(a,2)))
      integer::i
      do i=1,size(d);d(i)=a(i,i);end do
   end function diagonal

end module fhmm_estimation
