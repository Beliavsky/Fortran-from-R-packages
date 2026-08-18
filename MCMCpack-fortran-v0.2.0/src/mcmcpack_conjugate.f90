! SPDX-License-Identifier: GPL-3.0-only
! Direct posterior simulators translated from MCMCpack R/MCmodels.R.
module mcmcpack_conjugate
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : rbeta,rgamma_mt,rnorm,rdirichlet_rng
   implicit none
   private
   public :: mc_binomial_beta,mc_poisson_gamma,mc_normal_normal,mc_multinom_dirichlet
contains
   subroutine mc_binomial_beta(y,n,alpha,beta,draws,status)
      integer,intent(in)::y,n
      real(dp),intent(in)::alpha,beta
      real(dp),intent(out)::draws(:)
      integer,intent(out),optional::status
      integer::i
      if(present(status))status=0
      if(y<0.or.n<0.or.y>n.or.alpha<=0.0_dp.or.beta<=0.0_dp)then;if(present(status))status=1;draws=0.0_dp;return;end if
      do i=1,size(draws);draws(i)=rbeta(alpha+real(y,dp),beta+real(n-y,dp));end do
   end subroutine mc_binomial_beta

   subroutine mc_poisson_gamma(y,alpha,beta,draws,status)
      integer,intent(in)::y(:)
      real(dp),intent(in)::alpha,beta
      real(dp),intent(out)::draws(:)
      integer,intent(out),optional::status
      integer::i
      real(dp)::shape,rate
      if(present(status))status=0
      if(any(y<0).or.alpha<=0.0_dp.or.beta<=0.0_dp)then;if(present(status))status=1;draws=0.0_dp;return;end if
      shape=alpha+real(sum(y),dp);rate=beta+real(size(y),dp)
      do i=1,size(draws);draws(i)=rgamma_mt(shape,1.0_dp/rate);end do
   end subroutine mc_poisson_gamma

   subroutine mc_normal_normal(y,sigma2,mu0,tau20,draws,status)
      real(dp),intent(in)::y(:),sigma2,mu0,tau20
      real(dp),intent(out)::draws(:)
      integer,intent(out),optional::status
      integer::i
      real(dp)::mu1,tau21
      if(present(status))status=0
      if(sigma2<=0.0_dp.or.tau20<=0.0_dp)then;if(present(status))status=1;draws=0.0_dp;return;end if
      tau21=1.0_dp/(1.0_dp/tau20+real(size(y),dp)/sigma2)
      mu1=tau21*(mu0/tau20+sum(y)/sigma2)
      do i=1,size(draws);draws(i)=rnorm(mu1,sqrt(tau21));end do
   end subroutine mc_normal_normal

   subroutine mc_multinom_dirichlet(y,alpha0,draws,status)
      integer,intent(in)::y(:)
      real(dp),intent(in)::alpha0(:)
      real(dp),intent(out)::draws(:,:)
      integer,intent(out),optional::status
      real(dp)::a(size(alpha0))
      integer::i
      if(present(status))status=0
      if(size(y)/=size(alpha0).or.size(draws,2)/=size(alpha0).or.any(y<0).or.any(alpha0<=0.0_dp))then
         if(present(status))status=1;draws=0.0_dp;return
      end if
      a=alpha0+real(y,dp)
      do i=1,size(draws,1);call rdirichlet_rng(a,draws(i,:));end do
   end subroutine mc_multinom_dirichlet
end module mcmcpack_conjugate
