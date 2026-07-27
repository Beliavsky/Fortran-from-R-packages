! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_hierarchical
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use fhmm_kinds, only: dp, tiny_prob
   use fhmm_types, only: hhmm_parameters, hhmm_simulation, hhmm_inference_result
   use fhmm_math, only: stationary_distribution, log_sum_exp, seed_rng
   use fhmm_distributions, only: distribution_logpdf, distribution_random
   use fhmm_algorithms, only: hmm_log_likelihood, viterbi_decode, simulate_markov_chain
   use fhmm_parameters, only: validate_hhmm_parameters
   implicit none
   private
   public :: hhmm_log_likelihood, hhmm_forward_backward, simulate_hhmm_model
   public :: decode_hhmm

contains

   real(dp) function hhmm_log_likelihood(coarse_obs,fine_obs,lengths,par) result(ll)
      real(dp),intent(in)::coarse_obs(:),fine_obs(:,:)
      integer,intent(in)::lengths(:)
      type(hhmm_parameters),intent(in)::par
      type(hhmm_inference_result)::inf
      inf=hhmm_forward_backward(coarse_obs,fine_obs,lengths,par)
      if(inf%ok)then;ll=inf%log_likelihood;else;ll=-huge(1.0_dp);end if
   end function hhmm_log_likelihood

   function hhmm_forward_backward(coarse_obs,fine_obs,lengths,par) result(out)
      real(dp),intent(in)::coarse_obs(:),fine_obs(:,:)
      integer,intent(in)::lengths(:)
      type(hhmm_parameters),intent(in)::par
      type(hhmm_inference_result)::out
      real(dp),allocatable::logemit(:,:),phi(:,:),alpha(:,:),beta(:,:),delta(:),tmp(:)
      real(dp)::c,scale
      integer::m,Tn,t,s,j
      if(.not.validate_hhmm_parameters(par))then;out%message='invalid HHMM parameters';return;end if
      Tn=size(coarse_obs);m=size(par%coarse%mu)
      if(size(fine_obs,1)/=Tn.or.size(lengths)/=Tn)then;out%message='dimension mismatch';return;end if
      allocate(logemit(m,Tn),out%fine_log_likelihood(m,Tn),phi(m,Tn),tmp(m),delta(m))
      do t=1,Tn
         if(lengths(t)<1.or.lengths(t)>size(fine_obs,2))then;out%message='invalid chunk length';return;end if
         do s=1,m
            out%fine_log_likelihood(s,t)=hmm_log_likelihood(fine_obs(t,1:lengths(t)),par%fine(s))
            logemit(s,t)=distribution_logpdf(par%coarse%distribution,coarse_obs(t), &
               par%coarse%mu(s),par%coarse%sigma(s),par%coarse%df(s))+out%fine_log_likelihood(s,t)
         end do
      end do
      delta=stationary_distribution(par%coarse%gamma);phi(:,1)=log(max(delta,tiny_prob))+logemit(:,1)
      do t=2,Tn
         do j=1,m
            tmp=phi(:,t-1)+log(max(par%coarse%gamma(:,j),tiny_prob))
            phi(j,t)=log_sum_exp(tmp)+logemit(j,t)
         end do
      end do
      out%log_likelihood=log_sum_exp(phi(:,Tn))
      allocate(alpha(m,Tn),beta(m,Tn),out%coarse_filtered(m,Tn),out%coarse_smoothed(m,Tn))
      alpha(:,1)=delta*exp(logemit(:,1)-maxval(logemit(:,1)));scale=sum(alpha(:,1));alpha(:,1)=alpha(:,1)/scale
      out%coarse_filtered(:,1)=alpha(:,1)
      do t=2,Tn
         alpha(:,t)=matmul(transpose(par%coarse%gamma),alpha(:,t-1))* &
            exp(logemit(:,t)-maxval(logemit(:,t)))
         scale=sum(alpha(:,t));if(scale<=tiny_prob)then;out%message='zero forward scale';return;end if
         alpha(:,t)=alpha(:,t)/scale;out%coarse_filtered(:,t)=alpha(:,t)
      end do
      beta(:,Tn)=1.0_dp
      do t=Tn-1,1,-1
         tmp=exp(logemit(:,t+1)-maxval(logemit(:,t+1)))*beta(:,t+1)
         beta(:,t)=matmul(par%coarse%gamma,tmp)
         c=sum(beta(:,t));if(c>tiny_prob)beta(:,t)=beta(:,t)/c
      end do
      do t=1,Tn
         out%coarse_smoothed(:,t)=alpha(:,t)*beta(:,t);c=sum(out%coarse_smoothed(:,t))
         if(c>tiny_prob)out%coarse_smoothed(:,t)=out%coarse_smoothed(:,t)/c
      end do
      out%ok=.true.;out%message='ok'
   end function hhmm_forward_backward

   function simulate_hhmm_model(par,chunk_lengths,seed) result(out)
      type(hhmm_parameters),intent(in)::par
      integer,intent(in)::chunk_lengths(:)
      integer,intent(in),optional::seed
      type(hhmm_simulation)::out
      integer::Tn,maxlen,t,j,s
      real(dp)::nan
      if(present(seed))call seed_rng(seed)
      Tn=size(chunk_lengths);maxlen=maxval(chunk_lengths);nan=ieee_value(0.0_dp,ieee_quiet_nan)
      allocate(out%chunk_lengths(Tn),out%coarse_states(Tn))
      out%chunk_lengths=chunk_lengths
      out%coarse_states=simulate_markov_chain(par%coarse%gamma,Tn)
      allocate(out%coarse_observations(Tn),out%fine_states(Tn,maxlen),out%fine_observations(Tn,maxlen))
      out%fine_states=0;out%fine_observations=nan
      do t=1,Tn
         s=out%coarse_states(t)
         out%coarse_observations(t)=distribution_random(par%coarse%distribution,par%coarse%mu(s), &
            par%coarse%sigma(s),par%coarse%df(s))
         out%fine_states(t,1:chunk_lengths(t))=simulate_markov_chain(par%fine(s)%gamma,chunk_lengths(t))
         do j=1,chunk_lengths(t)
            out%fine_observations(t,j)=distribution_random(par%fine(s)%distribution, &
               par%fine(s)%mu(out%fine_states(t,j)),par%fine(s)%sigma(out%fine_states(t,j)), &
               par%fine(s)%df(out%fine_states(t,j)))
         end do
      end do
   end function simulate_hhmm_model

   subroutine decode_hhmm(coarse_obs,fine_obs,lengths,par,coarse_states,fine_states)
      real(dp),intent(in)::coarse_obs(:),fine_obs(:,:)
      integer,intent(in)::lengths(:)
      type(hhmm_parameters),intent(in)::par
      integer,allocatable,intent(out)::coarse_states(:),fine_states(:,:)
      type(hhmm_inference_result)::inf
            integer::t,s,m,maxlen
            ! Coarse Viterbi emissions include the fine-scale chunk likelihood.
      inf=hhmm_forward_backward(coarse_obs,fine_obs,lengths,par)
      m=size(par%coarse%mu);maxlen=size(fine_obs,2)
      allocate(coarse_states(size(coarse_obs)),fine_states(size(coarse_obs),maxlen));fine_states=0
      call coarse_viterbi(coarse_obs,inf%fine_log_likelihood,par,coarse_states)
      do t=1,size(coarse_obs)
         s=coarse_states(t)
         fine_states(t,1:lengths(t))=viterbi_decode(fine_obs(t,1:lengths(t)),par%fine(s))
      end do
   end subroutine decode_hhmm

   subroutine coarse_viterbi(obs,fine_ll,par,states)
      real(dp),intent(in)::obs(:),fine_ll(:,:)
      type(hhmm_parameters),intent(in)::par
      integer,intent(out)::states(:)
      real(dp),allocatable::xi(:,:),delta(:),tmp(:)
      integer,allocatable::back(:,:)
      integer::m,Tn,t,j,s
      m=size(par%coarse%mu);Tn=size(obs);allocate(xi(m,Tn),back(m,Tn),tmp(m));delta=stationary_distribution(par%coarse%gamma)
      do s=1,m
         xi(s,1)=log(max(delta(s),tiny_prob))+distribution_logpdf(par%coarse%distribution,obs(1), &
            par%coarse%mu(s),par%coarse%sigma(s),par%coarse%df(s))+fine_ll(s,1)
      end do
      do t=2,Tn
         do j=1,m
            tmp=xi(:,t-1)+log(max(par%coarse%gamma(:,j),tiny_prob));back(j,t)=maxloc(tmp,dim=1)
            xi(j,t)=maxval(tmp)+distribution_logpdf(par%coarse%distribution,obs(t),par%coarse%mu(j), &
               par%coarse%sigma(j),par%coarse%df(j))+fine_ll(j,t)
         end do
      end do
      states(Tn)=maxloc(xi(:,Tn),dim=1)
      do t=Tn-1,1,-1;states(t)=back(states(t+1),t+1);end do
   end subroutine coarse_viterbi

end module fhmm_hierarchical
