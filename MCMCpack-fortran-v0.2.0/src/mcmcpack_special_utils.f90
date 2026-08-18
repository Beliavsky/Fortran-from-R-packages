! SPDX-License-Identifier: GPL-3.0-only
! Shared numerical helpers for the specialized MCMCpack samplers.
module mcmcpack_special_utils
   use mcmcpack_kinds, only : dp, pi
   use mcmcpack_rng, only : runif, rdirichlet_rng
   implicit none
   private
   public :: sample_categorical, general_hmm_ffbs, transition_counts
   public :: draw_dirichlet_rows, normal_logpdf, poisson_logpmf
   public :: logistic_safe, logsumexp, normalize_prob
contains
   pure real(dp) function logistic_safe(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp+exp(-min(x,700.0_dp)))
      else
         p = exp(max(x,-700.0_dp))/(1.0_dp+exp(max(x,-700.0_dp)))
      end if
   end function logistic_safe

   pure real(dp) function normal_logpdf(x,mu,var) result(v)
      real(dp), intent(in) :: x,mu,var
      if (var <= 0.0_dp) then
         v=-huge(1.0_dp)
      else
         v=-0.5_dp*(log(2.0_dp*pi*var)+(x-mu)*(x-mu)/var)
      end if
   end function normal_logpdf

   pure real(dp) function poisson_logpmf(y,eta) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: eta
      if (y < 0 .or. eta > 700.0_dp) then
         v=-huge(1.0_dp)
      else
         v=real(y,dp)*eta-exp(eta)-log_gamma(real(y+1,dp))
      end if
   end function poisson_logpmf

   pure real(dp) function logsumexp(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x)==0) then
         v=-huge(1.0_dp); return
      end if
      m=maxval(x)
      if (m <= -0.5_dp*huge(1.0_dp)) then
         v=m
      else
         v=m+log(sum(exp(x-m)))
      end if
   end function logsumexp

   subroutine normalize_prob(w)
      real(dp), intent(inout) :: w(:)
      real(dp) :: z
      w=max(w,0.0_dp); z=sum(w)
      if (z <= tiny(1.0_dp)) then
         w=1.0_dp/real(max(1,size(w)),dp)
      else
         w=w/z
      end if
   end subroutine normalize_prob

   integer function sample_categorical(prob) result(draw)
      real(dp), intent(in) :: prob(:)
      real(dp) :: u,c,z
      integer :: j
      z=sum(max(prob,0.0_dp))
      if (z <= tiny(1.0_dp)) then
         draw=1+int(runif()*real(size(prob),dp)); draw=min(draw,size(prob)); return
      end if
      u=runif()*z;c=0.0_dp;draw=size(prob)
      do j=1,size(prob)
         c=c+max(prob(j),0.0_dp)
         if (u<=c) then; draw=j; return; end if
      end do
   end function sample_categorical

   subroutine general_hmm_ffbs(log_emit,p,pi0,state,prob_state,status)
      real(dp), intent(in) :: log_emit(:,:),p(:,:),pi0(:)
      integer, intent(out) :: state(size(log_emit,1))
      real(dp), intent(out) :: prob_state(size(log_emit,1),size(log_emit,2))
      integer, intent(out) :: status
      integer :: n,k,t,j,next
      real(dp), allocatable :: filt(:,:),pred(:),w(:),back(:,:),tmp(:)
      real(dp) :: mx,z
      n=size(log_emit,1);k=size(log_emit,2);status=0
      if(n<1.or.k<1.or.any(shape(p)/=[k,k]).or.size(pi0)/=k)then;status=1;return;end if
      allocate(filt(n,k),pred(k),w(k),back(n,k),tmp(k))
      pred=max(pi0,0.0_dp);call normalize_prob(pred)
      do t=1,n
         if(t>1)then;pred=matmul(filt(t-1,:),p);call normalize_prob(pred);end if
         mx=maxval(log_emit(t,:));w=pred*exp(log_emit(t,:)-mx);z=sum(w)
         if(z<=tiny(1.0_dp))then;status=2;return;end if
         filt(t,:)=w/z
      end do
      state(n)=sample_categorical(filt(n,:))
      do t=n-1,1,-1
         next=state(t+1);w=filt(t,:)*p(:,next);call normalize_prob(w);state(t)=sample_categorical(w)
      end do
      back(n,:)=1.0_dp
      prob_state(n,:)=filt(n,:);call normalize_prob(prob_state(n,:))
      do t=n-1,1,-1
         do j=1,k
            tmp=p(j,:)*exp(log_emit(t+1,:)-maxval(log_emit(t+1,:)))*back(t+1,:)
            back(t,j)=sum(tmp)
         end do
         if(maxval(back(t,:))>0.0_dp)back(t,:)=back(t,:)/maxval(back(t,:))
         prob_state(t,:)=filt(t,:)*back(t,:);call normalize_prob(prob_state(t,:))
      end do
   end subroutine general_hmm_ffbs

   subroutine transition_counts(state,k,count)
      integer, intent(in) :: state(:),k
      integer, intent(out) :: count(k,k)
      integer :: t
      count=0
      do t=1,size(state)-1
         if(state(t)>=1.and.state(t)<=k.and.state(t+1)>=1.and.state(t+1)<=k) &
            count(state(t),state(t+1))=count(state(t),state(t+1))+1
      end do
   end subroutine transition_counts

   subroutine draw_dirichlet_rows(alpha,p,status)
      real(dp), intent(in) :: alpha(:,:)
      real(dp), intent(out) :: p(size(alpha,1),size(alpha,2))
      integer, intent(out) :: status
      integer :: i,k
      k=size(alpha,2);status=0
      if(any(alpha<=0.0_dp))then;status=1;return;end if
      do i=1,size(alpha,1)
         call rdirichlet_rng(alpha(i,:),p(i,:))
      end do
   end subroutine draw_dirichlet_rows
end module mcmcpack_special_utils
