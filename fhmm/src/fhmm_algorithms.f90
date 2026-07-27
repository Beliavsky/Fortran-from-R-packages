! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_algorithms
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fhmm_kinds, only: dp, tiny_prob
   use fhmm_types, only: hmm_parameters, hmm_simulation, inference_result
   use fhmm_math, only: stationary_distribution, log_sum_exp, seed_rng
   use fhmm_distributions, only: distribution_logpdf, distribution_random
   use fhmm_parameters, only: validate_hmm_parameters
   implicit none
   private
   public :: emission_log_probabilities, hmm_log_likelihood, forward_backward
   public :: viterbi_decode, simulate_hmm_model, simulate_markov_chain

contains

   function emission_log_probabilities(observations,par) result(logp)
      real(dp),intent(in)::observations(:)
      type(hmm_parameters),intent(in)::par
      real(dp),allocatable::logp(:,:)
      integer::n,t,s
      n=size(par%mu);allocate(logp(n,size(observations)))
      do t=1,size(observations)
         do s=1,n
            logp(s,t)=distribution_logpdf(par%distribution,observations(t), &
               par%mu(s),par%sigma(s),par%df(s))
         end do
      end do
   end function emission_log_probabilities

   real(dp) function hmm_log_likelihood(observations,par) result(ll)
      real(dp),intent(in)::observations(:)
      type(hmm_parameters),intent(in)::par
      real(dp),allocatable::logp(:,:),phi(:,:),tmp(:),delta(:)
      integer::n,t,j
      if(.not.validate_hmm_parameters(par).or.size(observations)<1)then
         ll=-huge(1.0_dp);return
      end if
      n=size(par%mu);logp=emission_log_probabilities(observations,par)
      allocate(phi(n,size(observations)),tmp(n));delta=stationary_distribution(par%gamma)
      phi(:,1)=log(max(delta,tiny_prob))+logp(:,1)
      do t=2,size(observations)
         do j=1,n
            tmp=phi(:,t-1)+log(max(par%gamma(:,j),tiny_prob))
            phi(j,t)=log_sum_exp(tmp)+logp(j,t)
         end do
      end do
      tmp=phi(:,size(observations))
      ll=log_sum_exp(tmp)
   end function hmm_log_likelihood

   function forward_backward(observations,par) result(out)
      real(dp),intent(in)::observations(:)
      type(hmm_parameters),intent(in)::par
      type(inference_result)::out
      real(dp),allocatable::logp(:,:),alpha(:,:),beta(:,:),delta(:),b(:)
      real(dp)::scale
      integer::n,t,Tn
      if(.not.validate_hmm_parameters(par))then
         out%message='invalid HMM parameters';return
      end if
      Tn=size(observations);n=size(par%mu)
      if(Tn<1)then;out%message='no observations';return;end if
      logp=emission_log_probabilities(observations,par)
      allocate(alpha(n,Tn),beta(n,Tn),out%filtered(n,Tn),out%predicted(n,Tn), &
         out%smoothed(n,Tn),out%scales(Tn),delta(n),b(n))
      delta=stationary_distribution(par%gamma)
      out%predicted(:,1)=delta
      b=exp(max(logp(:,1),log(tiny_prob)))
      alpha(:,1)=delta*b;scale=sum(alpha(:,1))
      if(scale<=tiny_prob)then;out%message='zero forward scale';return;end if
      alpha(:,1)=alpha(:,1)/scale;out%scales(1)=scale;out%filtered(:,1)=alpha(:,1)
      out%log_likelihood=log(scale)
      do t=2,Tn
         out%predicted(:,t)=matmul(transpose(par%gamma),alpha(:,t-1))
         b=exp(max(logp(:,t),log(tiny_prob)))
         alpha(:,t)=out%predicted(:,t)*b;scale=sum(alpha(:,t))
         if(scale<=tiny_prob.or..not.ieee_is_finite(scale))then
            out%message='zero or nonfinite forward scale';return
         end if
         alpha(:,t)=alpha(:,t)/scale;out%scales(t)=scale
         out%filtered(:,t)=alpha(:,t);out%log_likelihood=out%log_likelihood+log(scale)
      end do
      beta(:,Tn)=1.0_dp
      do t=Tn-1,1,-1
         b=exp(max(logp(:,t+1),log(tiny_prob)))*beta(:,t+1)
         beta(:,t)=matmul(par%gamma,b)/out%scales(t+1)
      end do
      do t=1,Tn
         out%smoothed(:,t)=alpha(:,t)*beta(:,t)
         scale=sum(out%smoothed(:,t))
         if(scale>tiny_prob)out%smoothed(:,t)=out%smoothed(:,t)/scale
      end do
      out%ok=.true.;out%message='ok'
   end function forward_backward

   function viterbi_decode(observations,par) result(states)
      real(dp),intent(in)::observations(:)
      type(hmm_parameters),intent(in)::par
      integer,allocatable::states(:)
      real(dp),allocatable::logp(:,:),xi(:,:),delta(:),tmp(:)
      integer,allocatable::back(:,:)
      integer::n,t,j,Tn
      Tn=size(observations);n=size(par%mu);allocate(states(Tn))
      if(Tn==0.or..not.validate_hmm_parameters(par))then
         states=0;return
      end if
      logp=emission_log_probabilities(observations,par);delta=stationary_distribution(par%gamma)
      allocate(xi(n,Tn),back(n,Tn),tmp(n));xi(:,1)=log(max(delta,tiny_prob))+logp(:,1);back(:,1)=0
      do t=2,Tn
         do j=1,n
            tmp=xi(:,t-1)+log(max(par%gamma(:,j),tiny_prob))
            back(j,t)=maxloc(tmp,dim=1);xi(j,t)=maxval(tmp)+logp(j,t)
         end do
      end do
      states(Tn)=maxloc(xi(:,Tn),dim=1)
      do t=Tn-1,1,-1
         states(t)=back(states(t+1),t+1)
      end do
   end function viterbi_decode

   function simulate_markov_chain(gamma,n_steps,initial_state) result(states)
      real(dp),intent(in)::gamma(:,:)
      integer,intent(in)::n_steps
      integer,intent(in),optional::initial_state
      integer,allocatable::states(:)
      real(dp),allocatable::delta(:)
      real(dp)::u,c
      integer::t,j,n
      n=size(gamma,1);allocate(states(n_steps));delta=stationary_distribution(gamma)
      if(n_steps<1)return
      if(present(initial_state))then
         states(1)=initial_state
      else
         call random_number(u);c=0.0_dp;states(1)=n
         do j=1,n;c=c+delta(j);if(u<=c)then;states(1)=j;exit;end if;end do
      end if
      do t=2,n_steps
         call random_number(u);c=0.0_dp;states(t)=n
         do j=1,n;c=c+gamma(states(t-1),j);if(u<=c)then;states(t)=j;exit;end if;end do
      end do
   end function simulate_markov_chain

   function simulate_hmm_model(par,n_steps,seed,initial_state) result(out)
      type(hmm_parameters),intent(in)::par
      integer,intent(in)::n_steps
      integer,intent(in),optional::seed,initial_state
      type(hmm_simulation)::out
      integer::t,s
      if(present(seed))call seed_rng(seed)
      allocate(out%states(n_steps),out%observations(n_steps))
      if(present(initial_state))then
         out%states=simulate_markov_chain(par%gamma,n_steps,initial_state)
      else
         out%states=simulate_markov_chain(par%gamma,n_steps)
      end if
      do t=1,n_steps
         s=out%states(t)
         out%observations(t)=distribution_random(par%distribution,par%mu(s),par%sigma(s),par%df(s))
      end do
   end function simulate_hmm_model

end module fhmm_algorithms
