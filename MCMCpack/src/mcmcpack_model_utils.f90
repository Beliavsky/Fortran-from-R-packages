! SPDX-License-Identifier: GPL-3.0-only
! Computational summaries translated from BayesFactors.R, btsutil.R,
! hidden.R, make.breaklist.R and SSVSquantregsummary.R.
module mcmcpack_model_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mcmcpack_kinds, only : dp
   implicit none
   private
   public :: post_prob_mod, make_breaklist, agreement_matrix, transition_prior
   public :: marginal_inclusion, top_models, model_prob_result
   public :: bayes_factor, bayes_factor_result

   type :: bayes_factor_result
      real(dp), allocatable :: factor(:,:)
      real(dp), allocatable :: log_factor(:,:)
      real(dp), allocatable :: log_marginal(:)
   end type bayes_factor_result

   type :: model_prob_result
      integer, allocatable :: models(:,:)   ! unique 0/1 inclusion rows
      real(dp), allocatable :: probability(:)
   end type model_prob_result
contains
   function bayes_factor(log_marginal) result(out)
      ! Computational core of BayesFactor(): compare models fit to the same data.
      real(dp), intent(in) :: log_marginal(:)
      type(bayes_factor_result) :: out
      integer :: i,j,n
      n=size(log_marginal)
      allocate(out%factor(n,n),out%log_factor(n,n),out%log_marginal(n))
      out%log_marginal=log_marginal
      do i=1,n
         do j=1,n
            out%log_factor(i,j)=log_marginal(i)-log_marginal(j)
            if(out%log_factor(i,j)>log(huge(1.0_dp)))then
               out%factor(i,j)=huge(1.0_dp)
            else if(out%log_factor(i,j)<log(tiny(1.0_dp)))then
               out%factor(i,j)=0.0_dp
            else
               out%factor(i,j)=exp(out%log_factor(i,j))
            end if
         end do
      end do
   end function bayes_factor

   function post_prob_mod(log_marg,prior) result(prob)
      real(dp),intent(in)::log_marg(:)
      real(dp),intent(in),optional::prior(:)
      real(dp)::prob(size(log_marg)),lp(size(log_marg)),p(size(log_marg)),m,z
      integer::n
      n=size(log_marg);if(n==0)return
      if(present(prior))then
         if(size(prior)==1)then;p=prior(1)
         else if(size(prior)==n)then;p=prior
         else;prob=-1.0_dp;return;end if
      else;p=1.0_dp;end if
      if(any(p<=0.0_dp))then;prob=-1.0_dp;return;end if
      p=p/sum(p);lp=log_marg+log(p);m=maxval(lp);z=sum(exp(lp-m));prob=exp(lp-m)/z
   end function post_prob_mod

   function make_breaklist(logbf,threshold) result(out)
      real(dp),intent(in)::logbf(:,:)
      real(dp),intent(in),optional::threshold
      integer::out(size(logbf,1)),i,j,j1,j2
      real(dp)::th,best,second,ratio
      th=3.0_dp;if(present(threshold))th=threshold
      do i=1,size(logbf,1)
         if(any(.not.ieee_is_finite(logbf(i,:))))then;out(i)=0;cycle;end if
         j1=1;best=logbf(i,1)
         do j=2,size(logbf,2);if(logbf(i,j)>best)then;best=logbf(i,j);j1=j;end if;end do
         if(size(logbf,2)==1)then;out(i)=j1-1;cycle;end if
         j2=merge(2,1,j1==1);second=logbf(i,j2)
         do j=1,size(logbf,2);if(j/=j1.and.logbf(i,j)>second)then;second=logbf(i,j);j2=j;end if;end do
         ratio=exp(min(log(huge(1.0_dp)),best-second))
         if(ratio>th)then;out(i)=j1-1;else;out(i)=min(j1,j2)-1;end if
      end do
   end function make_breaklist

   function agreement_matrix(x,missing) result(a)
      integer,intent(in)::x(:,:)
      integer,intent(in),optional::missing
      real(dp)::a(size(x,1),size(x,1))
      integer::i,j,k,nitem,count
      nitem=size(x,2)
      if(present(missing)) continue
      do i=1,size(x,1);do j=1,size(x,1)
         count=0
         do k=1,nitem
            ! MCMCpack treats missingness as a category for agreement.
            if(x(i,k)==x(j,k))count=count+1
         end do
         a(i,j)=real(count,dp)/real(max(1,nitem),dp)
      end do;end do
   end function agreement_matrix

   function transition_prior(nbreak,nobs,a,b) result(trans)
      integer,intent(in)::nbreak,nobs
      real(dp),intent(in),optional::a,b
      real(dp)::trans(nbreak+1,nbreak+1),aa,bb
      integer::i
      bb=0.1_dp;if(present(b))bb=b
      aa=bb*real(nint(real(nobs,dp)/real(nbreak+1,dp)),dp);if(present(a))aa=a
      trans=0.0_dp
      do i=1,nbreak+1;trans(i,i)=1.0_dp;end do
      do i=1,nbreak;trans(i,i)=aa;trans(i,i+1)=bb;end do
   end function transition_prior

   function marginal_inclusion(gamma) result(prob)
      integer,intent(in)::gamma(:,:)
      real(dp)::prob(size(gamma,2))
      if(size(gamma,1)==0)then;prob=0.0_dp;else;prob=sum(real(gamma,dp),dim=1)/real(size(gamma,1),dp);end if
   end function marginal_inclusion

   function top_models(gamma,nmodels) result(out)
      integer,intent(in)::gamma(:,:)
      integer,intent(in),optional::nmodels
      type(model_prob_result)::out
      integer::n,k,nuniq,i,j,u,keep,im,jm
      integer,allocatable::uniq(:,:),counts(:),order(:)
      n=size(gamma,1);k=size(gamma,2);if(n==0)then;allocate(out%models(0,k),out%probability(0));return;end if
      allocate(uniq(n,k),counts(n));counts=0;nuniq=0
      do i=1,n
         u=0
         do j=1,nuniq;if(all(gamma(i,:)==uniq(j,:)))then;u=j;exit;end if;end do
         if(u==0)then;nuniq=nuniq+1;uniq(nuniq,:)=gamma(i,:);counts(nuniq)=1;else;counts(u)=counts(u)+1;end if
      end do
      allocate(order(nuniq));order=[(i,i=1,nuniq)]
      do i=1,nuniq-1
         im=i;do j=i+1,nuniq;if(counts(order(j))>counts(order(im)))im=j;end do
         if(im/=i)then;jm=order(i);order(i)=order(im);order(im)=jm;end if
      end do
      keep=min(nuniq,5);if(present(nmodels))keep=min(nuniq,max(0,nmodels))
      allocate(out%models(keep,k),out%probability(keep))
      do i=1,keep;out%models(i,:)=uniq(order(i),:);out%probability(i)=real(counts(order(i)),dp)/real(n,dp);end do
   end function top_models
end module mcmcpack_model_utils
