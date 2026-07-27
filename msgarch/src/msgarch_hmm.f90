! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_hmm
   use msgarch_kinds, only : dp, pi
   use msgarch_types, only : hmm_fit_result
   use msgarch_models, only : stationary_distribution, normalize_transition
   implicit none
   private
   public :: fit_gaussian_hmm, fit_gaussian_mixture, viterbi_gaussian_hmm
contains
   pure function gaussian_pdf(x,mean,variance) result(value)
      real(dp),intent(in)::x,mean,variance
      real(dp)::value
      value=exp(-0.5_dp*(x-mean)**2/variance)/sqrt(2.0_dp*pi*variance)
      value=max(value,1.0e-300_dp)
   end function gaussian_pdf

   function fit_gaussian_hmm(y,k,max_iterations,tolerance,zero_mean) result(fit)
      real(dp),intent(in)::y(:)
      integer,intent(in)::k
      integer,intent(in),optional::max_iterations
      real(dp),intent(in),optional::tolerance
      logical,intent(in),optional::zero_mean
      type(hmm_fit_result)::fit
      real(dp),allocatable::alpha(:,:),beta(:,:),gamma(:,:),xi(:,:,:),emission(:,:),scale(:),delta(:)
      real(dp),allocatable::old_mean(:),old_variance(:),old_transition(:,:)
      real(dp)::overall_mean,overall_var,tol,denom,loglik_old,loglik
      integer::t,i,j,maxit,iter,n
      logical::force_zero
      n=size(y);maxit=1000;if(present(max_iterations))maxit=max_iterations
      tol=1.0e-8_dp;if(present(tolerance))tol=tolerance
      force_zero=.false.;if(present(zero_mean))force_zero=zero_mean
      allocate(fit%mean(k),fit%variance(k),fit%transition(k,k),fit%probability(k))
      allocate(alpha(n,k),beta(n,k),gamma(n,k),xi(n-1,k,k),emission(n,k),scale(n),delta(k))
      allocate(old_mean(k),old_variance(k),old_transition(k,k))
      overall_mean=sum(y)/real(n,dp);overall_var=sum((y-overall_mean)**2)/real(max(n-1,1),dp)
      do i=1,k
         if(force_zero)then;fit%mean(i)=0.0_dp;else;fit%mean(i)=overall_mean*(0.8_dp+0.4_dp*real(i-1,dp)/real(max(k-1,1),dp));end if
         fit%variance(i)=overall_var*(0.6_dp+0.8_dp*real(i-1,dp)/real(max(k-1,1),dp))
      end do
      if(k==1)then;fit%transition=1.0_dp;else
         fit%transition=0.1_dp/real(k-1,dp);do i=1,k;fit%transition(i,i)=0.9_dp;end do
      end if
      loglik_old=-huge(1.0_dp)
      do iter=1,maxit
         old_mean=fit%mean;old_variance=fit%variance;old_transition=fit%transition
         delta=stationary_distribution(fit%transition)
         do t=1,n;do i=1,k;emission(t,i)=gaussian_pdf(y(t),fit%mean(i),fit%variance(i));end do;end do
         alpha(1,:)=delta*emission(1,:);scale(1)=sum(alpha(1,:));alpha(1,:)=alpha(1,:)/scale(1)
         do t=2,n
            alpha(t,:)=matmul(alpha(t-1,:),fit%transition)*emission(t,:)
            scale(t)=sum(alpha(t,:));alpha(t,:)=alpha(t,:)/scale(t)
         end do
         loglik=sum(log(scale));beta(n,:)=1.0_dp
         do t=n-1,1,-1
            beta(t,:)=matmul(fit%transition,emission(t+1,:)*beta(t+1,:))/scale(t+1)
         end do
         gamma=alpha*beta
         do t=1,n;gamma(t,:)=gamma(t,:)/sum(gamma(t,:));end do
         do t=1,n-1
            denom=0.0_dp
            do i=1,k;do j=1,k
               xi(t,i,j)=alpha(t,i)*fit%transition(i,j)*emission(t+1,j)*beta(t+1,j);denom=denom+xi(t,i,j)
            end do;end do
            xi(t,:,:)=xi(t,:,:)/max(denom,tiny(1.0_dp))
         end do
         do i=1,k
            denom=sum(gamma(:,i))
            if(.not.force_zero)fit%mean(i)=sum(gamma(:,i)*y)/max(denom,tiny(1.0_dp))
            fit%variance(i)=sum(gamma(:,i)*(y-fit%mean(i))**2)/max(denom,tiny(1.0_dp))
            fit%variance(i)=max(fit%variance(i),1.0e-10_dp)
            if(n>1)then
               do j=1,k;fit%transition(i,j)=sum(xi(:,i,j))/max(sum(gamma(1:n-1,i)),tiny(1.0_dp));end do
            end if
         end do
         call normalize_transition(fit%transition)
         fit%iterations=iter;fit%loglik=loglik
         if(abs(loglik-loglik_old)<tol*(1.0_dp+abs(loglik)))then;fit%converged=.true.;exit;end if
         loglik_old=loglik
      end do
      fit%probability=stationary_distribution(fit%transition)
   end function fit_gaussian_hmm

   function fit_gaussian_mixture(y,k,max_iterations,tolerance,zero_mean) result(fit)
      real(dp),intent(in)::y(:)
      integer,intent(in)::k
      integer,intent(in),optional::max_iterations
      real(dp),intent(in),optional::tolerance
      logical,intent(in),optional::zero_mean
      type(hmm_fit_result)::fit
      real(dp),allocatable::weight(:,:),old_mean(:),old_variance(:),old_probability(:)
      real(dp)::overall_mean,overall_var,tol,total_density,loglik,old_loglik,denom
      integer::n,i,t,maxit,iter
      logical::force_zero
      n=size(y);maxit=1000;if(present(max_iterations))maxit=max_iterations
      tol=1.0e-8_dp;if(present(tolerance))tol=tolerance
      force_zero=.false.;if(present(zero_mean))force_zero=zero_mean
      allocate(fit%mean(k),fit%variance(k),fit%probability(k),fit%transition(k,k),weight(n,k))
      allocate(old_mean(k),old_variance(k),old_probability(k))
      overall_mean=sum(y)/real(n,dp);overall_var=sum((y-overall_mean)**2)/real(max(n-1,1),dp)
      fit%probability=1.0_dp/real(k,dp)
      do i=1,k
         if(force_zero)then;fit%mean(i)=0.0_dp;else;fit%mean(i)=overall_mean*(0.8_dp+0.4_dp*real(i-1,dp)/real(max(k-1,1),dp));end if
         fit%variance(i)=overall_var*(0.6_dp+0.8_dp*real(i-1,dp)/real(max(k-1,1),dp))
      end do
      old_loglik=-huge(1.0_dp)
      do iter=1,maxit
         old_mean=fit%mean;old_variance=fit%variance;old_probability=fit%probability;loglik=0.0_dp
         do t=1,n
            do i=1,k;weight(t,i)=fit%probability(i)*gaussian_pdf(y(t),fit%mean(i),fit%variance(i));end do
            total_density=sum(weight(t,:));weight(t,:)=weight(t,:)/total_density;loglik=loglik+log(total_density)
         end do
         do i=1,k
            denom=sum(weight(:,i));fit%probability(i)=denom/real(n,dp)
            if(.not.force_zero)fit%mean(i)=sum(weight(:,i)*y)/max(denom,tiny(1.0_dp))
            fit%variance(i)=sum(weight(:,i)*(y-fit%mean(i))**2)/max(denom,tiny(1.0_dp));fit%variance(i)=max(fit%variance(i),1.0e-10_dp)
         end do
         fit%iterations=iter;fit%loglik=loglik
         if(abs(loglik-old_loglik)<tol*(1.0_dp+abs(loglik)))then;fit%converged=.true.;exit;end if
         old_loglik=loglik
      end do
      do i=1,k;fit%transition(i,:)=fit%probability;end do
   end function fit_gaussian_mixture

   function viterbi_gaussian_hmm(y,mean,variance,transition) result(states)
      real(dp),intent(in)::y(:),mean(:),variance(:),transition(:,:)
      integer,allocatable::states(:),back(:,:)
      real(dp),allocatable::score(:,:),delta(:)
      real(dp)::candidate,best
      integer::n,k,t,i,j
      n=size(y);k=size(mean);allocate(states(n),back(n,k),score(n,k));delta=stationary_distribution(transition)
      do i=1,k;score(1,i)=log(max(delta(i),tiny(1.0_dp)))+log(gaussian_pdf(y(1),mean(i),variance(i)));back(1,i)=1;end do
      do t=2,n
         do j=1,k
            best=-huge(1.0_dp);back(t,j)=1
            do i=1,k
               candidate=score(t-1,i)+log(max(transition(i,j),tiny(1.0_dp)))
               if(candidate>best)then;best=candidate;back(t,j)=i;end if
            end do
            score(t,j)=best+log(gaussian_pdf(y(t),mean(j),variance(j)))
         end do
      end do
      states(n)=maxloc(score(n,:),dim=1)
      do t=n,2,-1;states(t-1)=back(t,states(t));end do
   end function viterbi_gaussian_hmm
end module msgarch_hmm
