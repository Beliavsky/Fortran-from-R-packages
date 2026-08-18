! SPDX-License-Identifier: GPL-3.0-only
! Two-dimensional paired-comparison sampler translated from cMCMCpaircompare2d.cc.
module mcmcpack_paircompare2d
   use mcmcpack_kinds, only : dp, pi
   use mcmcpack_rng, only : runif, rnorm, rtruncnorm, rmvnorm
   use mcmcpack_linalg, only : inv_spd
   implicit none
   private

   type, public :: paircompare2d_result
      real(dp), allocatable :: draws(:,:)
      real(dp), allocatable :: gamma_accept_rate(:)
      integer :: status = 0
   end type paircompare2d_result

   public :: mcmc_paircompare2d, update_theta_candidate
contains
   function mcmc_paircompare2d(md,theta_start,gamma_start,theta_eq,theta_ineq,tune, &
                               burnin,mcmc,thin,store_theta,store_gamma) result(res)
      ! md(:,1) = respondent/judge, md(:,2) = candidate 1,
      ! md(:,3) = candidate 2, md(:,4) = chosen candidate.  IDs are 1-based.
      integer, intent(in) :: md(:,:)
      real(dp), intent(in) :: theta_start(:,:), gamma_start(:)
      real(dp), intent(in) :: theta_eq(:,:), theta_ineq(:,:), tune
      integer, intent(in) :: burnin,mcmc,thin
      logical, intent(in), optional :: store_theta,store_gamma
      type(paircompare2d_result) :: res
      integer :: n,ncand,nresp,nsamp,ncol,iter,count,i,j,r,c1,c2,info,col
      real(dp) :: theta(size(theta_start,1),2),gamma(size(gamma_start))
      real(dp), allocatable :: ystar(:),gamma_trial(:),gamma_accept(:)
      real(dp) :: mu,gold,gnew,llold,llnew,u,eta,etanew
      logical :: st,sg

      st=.true.; sg=.true.
      if(present(store_theta)) st=store_theta
      if(present(store_gamma)) sg=store_gamma
      n=size(md,1); ncand=size(theta_start,1); nresp=size(gamma_start)
      if(thin<=0.or.mcmc<=0)then;res%status=1;return;end if
      nsamp=mcmc/thin
      if(size(md,2)<4.or.size(theta_start,2)/=2.or.any(shape(theta_eq)/=[ncand,2]).or. &
         any(shape(theta_ineq)/=[ncand,2]).or.nsamp<=0.or.(.not.st.and..not.sg).or.tune<=0.0_dp)then
         res%status=1;return
      end if
      if(any(md(:,1)<1).or.any(md(:,1)>nresp).or.any(md(:,2)<1).or.any(md(:,2)>ncand).or. &
         any(md(:,3)<1).or.any(md(:,3)>ncand))then;res%status=2;return;end if
      if(any(gamma_start<0.0_dp).or.any(gamma_start>0.5_dp*pi))then;res%status=3;return;end if

      ncol=merge(2*ncand,0,st)+merge(nresp,0,sg)
      allocate(res%draws(nsamp,ncol),res%gamma_accept_rate(nresp),ystar(n), &
               gamma_trial(nresp),gamma_accept(nresp))
      theta=theta_start;gamma=gamma_start;gamma_trial=0.0_dp;gamma_accept=0.0_dp;count=0

      do iter=0,burnin+mcmc-1
         ! Latent Gaussian utilities, conditional on the observed winner.
         do i=1,n
            r=md(i,1);c1=md(i,2);c2=md(i,3)
            mu=cos(gamma(r))*(theta(c1,1)-theta(c2,1))+ &
               sin(gamma(r))*(theta(c1,2)-theta(c2,2))
            if(md(i,4)==c1)then
               ystar(i)=rtruncnorm(mu,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else if(md(i,4)==c2)then
               ystar(i)=rtruncnorm(mu,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
            else
               ystar(i)=rnorm(mu,1.0_dp)
            end if
         end do

         ! Respondent angles: uniform random-walk MH on [0,pi/2], with a
         ! uniform prior exactly as in the original implementation.
         do r=1,nresp
            gold=gamma(r)
            do
               gnew=gold+(1.0_dp-2.0_dp*runif())*tune
               if(gnew>=0.0_dp.and.gnew<=0.5_dp*pi)exit
            end do
            llold=0.0_dp;llnew=0.0_dp
            do i=1,n
               if(md(i,1)/=r)cycle
               c1=md(i,2);c2=md(i,3)
               eta=cos(gold)*(theta(c1,1)-theta(c2,1))+sin(gold)*(theta(c1,2)-theta(c2,2))
               etanew=cos(gnew)*(theta(c1,1)-theta(c2,1))+sin(gnew)*(theta(c1,2)-theta(c2,2))
               llold=llold-0.5_dp*(ystar(i)-eta)**2
               llnew=llnew-0.5_dp*(ystar(i)-etanew)**2
            end do
            gamma_trial(r)=gamma_trial(r)+1.0_dp
            u=runif()
            if(log(u)<min(0.0_dp,llnew-llold))then
               gamma(r)=gnew;gamma_accept(r)=gamma_accept(r)+1.0_dp
            end if
         end do

         ! Candidate positions, including MCMCpack's -999 equality sentinel
         ! and signed inequality constraints.
         do j=1,ncand
            call update_theta_candidate(j,md,ystar,gamma,theta,theta_eq,theta_ineq,info)
            if(info/=0)then;res%status=10+info;return;end if
         end do

         if(iter>=burnin.and.mod(iter,thin)==0)then
            count=count+1;col=0
            if(st)then
               res%draws(count,1:ncand)=theta(:,1)
               res%draws(count,ncand+1:2*ncand)=theta(:,2)
               col=2*ncand
            end if
            if(sg)res%draws(count,col+1:col+nresp)=gamma
         end if
      end do
      do r=1,nresp
         if(gamma_trial(r)>0.0_dp)then
            res%gamma_accept_rate(r)=gamma_accept(r)/gamma_trial(r)
         else
            res%gamma_accept_rate(r)=0.0_dp
         end if
      end do
   end function mcmc_paircompare2d

   subroutine update_theta_candidate(j,md,ystar,gamma,theta,theta_eq,theta_ineq,info)
      integer,intent(in)::j,md(:,:)
      real(dp),intent(in)::ystar(:),gamma(:),theta_eq(:,:),theta_ineq(:,:)
      real(dp),intent(inout)::theta(:,:)
      integer,intent(out)::info
      integer::i,r,c1,c2,nobs,ncand,k
      real(dp),allocatable::x(:,:),z(:)
      real(dp)::sgn,prec1,rhs1,v,m,lo,hi
      real(dp)::prec(2,2),cov(2,2),rhs(2),mean(2),draw(2)
      logical::free1,free2,ok

      info=0;ncand=size(theta,1);nobs=0
      do i=1,size(md,1)
         if(md(i,2)==j.or.md(i,3)==j)nobs=nobs+1
      end do
      allocate(x(nobs,2),z(nobs));k=0
      do i=1,size(md,1)
         c1=md(i,2);c2=md(i,3)
         if(c1/=j.and.c2/=j)cycle
         k=k+1;r=md(i,1)
         if(c1==j)then
            sgn=1.0_dp
            x(k,1)=cos(gamma(r));x(k,2)=sin(gamma(r))
            z(k)=ystar(i)+cos(gamma(r))*theta(c2,1)+sin(gamma(r))*theta(c2,2)
         else
            sgn=-1.0_dp
            x(k,1)=-cos(gamma(r));x(k,2)=-sin(gamma(r))
            z(k)=ystar(i)-cos(gamma(r))*theta(c1,1)-sin(gamma(r))*theta(c1,2)
         end if
         if(abs(sgn)>2.0_dp)info=99 ! unreachable; retains explicit sign semantics.
      end do

      free1=theta_eq(j,1)<-998.5_dp;free2=theta_eq(j,2)<-998.5_dp
      if(.not.free1.and..not.free2)then
         theta(j,:)=theta_eq(j,:);return
      end if
      if(.not.free1)theta(j,1)=theta_eq(j,1)
      if(.not.free2)theta(j,2)=theta_eq(j,2)

      if(free1.and.free2)then
         prec=matmul(transpose(x),x)
         prec(1,1)=prec(1,1)+1.0_dp;prec(2,2)=prec(2,2)+1.0_dp
         rhs=matmul(transpose(x),z)
         call inv_spd(prec,cov,info);if(info/=0)return
         mean=matmul(cov,rhs)
         do
            call rmvnorm(mean,cov,draw,info);if(info/=0)return
            ok=.true.
            if(theta_ineq(j,1)/=0.0_dp)ok=ok.and.theta_ineq(j,1)*draw(1)>=0.0_dp
            if(theta_ineq(j,2)/=0.0_dp)ok=ok.and.theta_ineq(j,2)*draw(2)>=0.0_dp
            if(ok)exit
         end do
         theta(j,:)=draw
      else if(free1)then
         prec1=1.0_dp;rhs1=0.0_dp
         do k=1,nobs
            prec1=prec1+x(k,1)*x(k,1)
            rhs1=rhs1+x(k,1)*(z(k)-x(k,2)*theta(j,2))
         end do
         v=1.0_dp/prec1;m=v*rhs1
         call constraint_bounds(theta_ineq(j,1),lo,hi)
         if(theta_ineq(j,1)==0.0_dp)then;theta(j,1)=rnorm(m,sqrt(v))
         else;theta(j,1)=rtruncnorm(m,sqrt(v),lo,hi);end if
      else
         prec1=1.0_dp;rhs1=0.0_dp
         do k=1,nobs
            prec1=prec1+x(k,2)*x(k,2)
            rhs1=rhs1+x(k,2)*(z(k)-x(k,1)*theta(j,1))
         end do
         v=1.0_dp/prec1;m=v*rhs1
         call constraint_bounds(theta_ineq(j,2),lo,hi)
         if(theta_ineq(j,2)==0.0_dp)then;theta(j,2)=rnorm(m,sqrt(v))
         else;theta(j,2)=rtruncnorm(m,sqrt(v),lo,hi);end if
      end if
      if(any(abs(theta(j,:))>huge(1.0_dp)/100.0_dp))info=98
      if(ncand<1)info=97
   end subroutine update_theta_candidate

   subroutine constraint_bounds(sign_constraint,lo,hi)
      real(dp),intent(in)::sign_constraint
      real(dp),intent(out)::lo,hi
      real(dp),parameter::big=huge(1.0_dp)/10.0_dp
      if(sign_constraint>0.0_dp)then;lo=0.0_dp;hi=big
      else if(sign_constraint<0.0_dp)then;lo=-big;hi=0.0_dp
      else;lo=-big;hi=big;end if
   end subroutine constraint_bounds
end module mcmcpack_paircompare2d
