! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_diagnostics
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use fhmm_kinds, only: dp, tiny_prob
   use fhmm_types, only: hmm_parameters, forecast_result, model_comparison
   use fhmm_math, only: normal_quantile
   use fhmm_distributions, only: distribution_cdf,distribution_quantile,distribution_mean
   use fhmm_algorithms, only: viterbi_decode
   implicit none
   private
   public :: pseudo_residuals, forecast_hmm, reorder_hmm_states, compare_hmm_model

contains

   function pseudo_residuals(observations,states,par) result(residuals)
      real(dp),intent(in)::observations(:)
      integer,intent(in)::states(:)
      type(hmm_parameters),intent(in)::par
      real(dp),allocatable::residuals(:)
      real(dp)::p
      integer::t,s
      allocate(residuals(size(observations)))
      if(size(states)/=size(observations))then
         residuals=ieee_value(0.0_dp,ieee_quiet_nan);return
      end if
      do t=1,size(observations)
         s=states(t);p=distribution_cdf(par%distribution,observations(t),par%mu(s),par%sigma(s),par%df(s))
         p=min(max(p,1.0e-12_dp),1.0_dp-1.0e-12_dp);residuals(t)=normal_quantile(p)
      end do
   end function pseudo_residuals

   function forecast_hmm(par,ahead,alpha,last_state,last_probabilities,upstream_quantiles) result(out)
      type(hmm_parameters),intent(in)::par
      integer,intent(in)::ahead
      real(dp),intent(in),optional::alpha,last_probabilities(:)
      integer,intent(in),optional::last_state
      logical,intent(in),optional::upstream_quantiles
      type(forecast_result)::out
      real(dp),allocatable::prob(:)
      real(dp)::a
      integer::h,n
      logical::compat
      n=size(par%mu);a=0.05_dp;if(present(alpha))a=alpha;compat=.false.;if(present(upstream_quantiles))compat=upstream_quantiles
      if(ahead<1.or.a<=0.0_dp.or.a>=0.5_dp)then;out%message='invalid forecast arguments';return;end if
      allocate(prob(n),out%state_probabilities(n,ahead),out%lower(ahead),out%median(ahead),out%upper(ahead),out%mean(ahead))
      if(present(last_probabilities))then
         if(size(last_probabilities)/=n)then;out%message='probability dimension mismatch';return;end if
         prob=last_probabilities/sum(last_probabilities)
      else if(present(last_state))then
         prob=0.0_dp;prob(last_state)=1.0_dp
      else
         prob=1.0_dp/real(n,dp)
      end if
      do h=1,ahead
         prob=matmul(transpose(par%gamma),prob);out%state_probabilities(:,h)=prob
         out%mean(h)=sum(prob*[(distribution_mean(par%distribution,par%mu(n),par%sigma(n),par%df(n)),n=1,size(par%mu))])
         if(compat)then
            out%lower(h)=weighted_component_quantile(a,prob,par)
            out%median(h)=weighted_component_quantile(0.5_dp,prob,par)
            out%upper(h)=weighted_component_quantile(1.0_dp-a,prob,par)
         else
            out%lower(h)=mixture_quantile(a,prob,par)
            out%median(h)=mixture_quantile(0.5_dp,prob,par)
            out%upper(h)=mixture_quantile(1.0_dp-a,prob,par)
         end if
      end do
      out%ok=.true.;out%message='ok'
   end function forecast_hmm

   pure real(dp) function weighted_component_quantile(p,prob,par) result(q)
      real(dp),intent(in)::p,prob(:)
      type(hmm_parameters),intent(in)::par
      integer::s
      q=0.0_dp
      do s=1,size(prob);q=q+prob(s)*distribution_quantile(par%distribution,p,par%mu(s),par%sigma(s),par%df(s));end do
   end function weighted_component_quantile

   pure real(dp) function mixture_quantile(p,prob,par) result(q)
      real(dp),intent(in)::p,prob(:)
      type(hmm_parameters),intent(in)::par
      real(dp)::lo,hi,mid,c
      integer::s,it
      lo=huge(1.0_dp);hi=-huge(1.0_dp)
      do s=1,size(prob)
         lo=min(lo,distribution_quantile(par%distribution,1.0e-8_dp,par%mu(s),par%sigma(s),par%df(s)))
         hi=max(hi,distribution_quantile(par%distribution,1.0_dp-1.0e-8_dp,par%mu(s),par%sigma(s),par%df(s)))
      end do
      do it=1,200
         mid=0.5_dp*(lo+hi);c=0.0_dp
         do s=1,size(prob);c=c+prob(s)*distribution_cdf(par%distribution,mid,par%mu(s),par%sigma(s),par%df(s));end do
         if(c<p)then;lo=mid;else;hi=mid;end if
      end do
      q=0.5_dp*(lo+hi)
   end function mixture_quantile

   function reorder_hmm_states(par,order) result(out)
      type(hmm_parameters),intent(in)::par
      integer,intent(in),optional::order(:)
      type(hmm_parameters)::out
      integer,allocatable::ord(:)
      integer::n,i,j
      n=size(par%mu);allocate(ord(n))
      if(present(order))then;ord=order;else;ord=[(i,i=1,n)];call sort_order(par%mu,ord);end if
      out%distribution=par%distribution;allocate(out%gamma(n,n),out%mu(n),out%sigma(n),out%df(n))
      do i=1,n
         out%mu(i)=par%mu(ord(i));out%sigma(i)=par%sigma(ord(i));out%df(i)=par%df(ord(i))
         do j=1,n;out%gamma(i,j)=par%gamma(ord(i),ord(j));end do
      end do
   end function reorder_hmm_states

   pure subroutine sort_order(x,ord)
      real(dp),intent(in)::x(:)
      integer,intent(inout)::ord(:)
      integer::i,j,v
      do i=2,size(ord);v=ord(i);j=i-1;do while(j>=1);if(x(ord(j))<=x(v))exit;ord(j+1)=ord(j);j=j-1;end do;ord(j+1)=v;end do
   end subroutine sort_order

   pure function compare_hmm_model(log_likelihood,n_parameters,n_observations) result(out)
      real(dp),intent(in)::log_likelihood
      integer,intent(in)::n_parameters,n_observations
      type(model_comparison)::out
      out%log_likelihood=log_likelihood;out%parameters=n_parameters;out%observations=n_observations
      out%aic=-2.0_dp*log_likelihood+2.0_dp*real(n_parameters,dp)
      out%bic=-2.0_dp*log_likelihood+real(n_parameters,dp)*log(real(n_observations,dp))
   end function compare_hmm_model

end module fhmm_diagnostics
