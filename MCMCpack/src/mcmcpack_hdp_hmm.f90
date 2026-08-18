! SPDX-License-Identifier: GPL-3.0-only
! Finite-truncation sticky HDP-HMM / HDP-HSMM count samplers.
! Based on cHDPHMMpoisson.cc, cHDPHMMnegbin.cc and cHDPHSMMnegbin.cc.
! The finite K weak-limit representation and CRT auxiliary-variable updates are explicit.
module mcmcpack_hdp_hmm
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : runif, rnorm, rbeta, rgamma_mt, rdirichlet_rng
   use mcmcpack_math, only : logsumexp
   use mcmcpack_special_utils, only : general_hmm_ffbs, poisson_logpmf
   implicit none
   private
   public :: hdp_count_result, hdphmm_poisson, hdphmm_negbin, hdphsmm_negbin

   type :: hdp_count_result
      real(dp),allocatable :: beta(:,:,:)
      real(dp),allocatable :: p(:,:,:)
      integer,allocatable :: state(:,:)
      real(dp),allocatable :: state_prob(:,:,:)
      real(dp),allocatable :: rho(:,:)
      real(dp),allocatable :: omega(:,:)
      real(dp),allocatable :: gamma(:),alpha_kappa(:),theta(:)
      integer :: status=0
   end type hdp_count_result
contains
   pure real(dp) function quad_prior(b,b0,bprec) result(v)
      real(dp),intent(in)::b(:),b0(:),bprec(:,:)
      real(dp)::d(size(b))
      d=b-b0;v=-0.5_dp*dot_product(d,matmul(bprec,d))
   end function quad_prior

   pure real(dp) function nb_logpmf(y,mu,rho) result(v)
      integer,intent(in)::y
      real(dp),intent(in)::mu,rho
      if(y<0.or.mu<=0.0_dp.or.rho<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      v=log_gamma(rho+real(y,dp))-log_gamma(rho)-log_gamma(real(y+1,dp))+rho*log(rho)+ &
        real(y,dp)*log(mu)-(rho+real(y,dp))*log(rho+mu)
   end function nb_logpmf

   pure real(dp) function state_beta_logpost(y,x,state,j,b,b0,bprec,family,rho) result(v)
      integer,intent(in)::y(:),state(:),j,family
      real(dp),intent(in)::x(:,:),b(:),b0(:),bprec(:,:),rho
      integer::t
      real(dp)::mu
      v=quad_prior(b,b0,bprec)
      do t=1,size(y)
         if(state(t)/=j)cycle;mu=exp(max(-700.0_dp,min(700.0_dp,dot_product(x(t,:),b))))
         if(family==1)then;v=v+poisson_logpmf(y(t),mu)
         else;v=v+nb_logpmf(y(t),mu,rho);end if
      end do
   end function state_beta_logpost

   subroutine beta_rw_update(y,x,state,j,b,b0,bprec,family,rho,step)
      integer,intent(in)::y(:),state(:),j,family
      real(dp),intent(in)::x(:,:),b0(:),bprec(:,:),rho,step
      real(dp),intent(inout)::b(:)
      integer::h
      real(dp)::cand(size(b)),lo,ln
      do h=1,size(b)
         cand=b;cand(h)=b(h)+step*rnorm()
         lo=state_beta_logpost(y,x,state,j,b,b0,bprec,family,rho)
         ln=state_beta_logpost(y,x,state,j,cand,b0,bprec,family,rho)
         if(log(max(runif(),tiny(1.0_dp)))<min(0.0_dp,ln-lo))b=cand
      end do
   end subroutine beta_rw_update

   pure real(dp) function rho_logpost(y,x,state,j,b,rho,g,e,f) result(v)
      integer,intent(in)::y(:),state(:),j
      real(dp),intent(in)::x(:,:),b(:),rho,g,e,f
      integer::t
      real(dp)::mu
      if(rho<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      v=(e-1.0_dp)*log(rho)-(e+f)*log(rho+g)
      do t=1,size(y)
         if(state(t)/=j)cycle;mu=exp(max(-700.0_dp,min(700.0_dp,dot_product(x(t,:),b))))
         v=v+nb_logpmf(y(t),mu,rho)
      end do
   end function rho_logpost

   subroutine rho_slice_update(y,x,state,j,b,rho,step,g,e,f)
      integer,intent(in)::y(:),state(:),j
      real(dp),intent(in)::x(:,:),b(:),step,g,e,f
      real(dp),intent(inout)::rho
      integer::jl,jr,it
      real(dp)::slice,left,right,cand,lc
      slice=rho_logpost(y,x,state,j,b,rho,g,e,f)+log(max(runif(),tiny(1.0_dp)))
      left=max(1.0e-10_dp,rho-step*runif());right=left+step;jl=int(100*runif());jr=99-jl
      do while(jl>0.and.rho_logpost(y,x,state,j,b,left,g,e,f)>slice)
         left=max(1.0e-10_dp,left-step);jl=jl-1;if(left<=1.0e-10_dp)exit
      end do
      do while(jr>0.and.rho_logpost(y,x,state,j,b,right,g,e,f)>slice);right=right+step;jr=jr-1;end do
      do it=1,10000
         cand=left+(right-left)*runif();lc=rho_logpost(y,x,state,j,b,cand,g,e,f)
         if(lc>slice)then;rho=cand;return;end if
         if(cand<rho)then;left=cand;else;right=cand;end if
      end do
   end subroutine rho_slice_update

   integer function crt_draw(n,a) result(m)
      integer,intent(in)::n
      real(dp),intent(in)::a
      integer::i
      m=0;if(n<=0.or.a<=0.0_dp)return
      do i=1,n;if(runif()<a/(a+real(i-1,dp)))m=m+1;end do
   end function crt_draw

   subroutine transition_counts(state,k,ntr)
      integer,intent(in)::state(:),k
      integer,intent(out)::ntr(k,k)
      integer::t
      ntr=0
      do t=1,size(state)-1
         if(state(t)>=1.and.state(t)<=k.and.state(t+1)>=1.and.state(t+1)<=k) &
            ntr(state(t),state(t+1))=ntr(state(t),state(t+1))+1
      end do
   end subroutine transition_counts

   subroutine sticky_hdp_update(state,p,gw,gamma,ak,theta,a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta)
      integer,intent(in)::state(:)
      real(dp),intent(inout)::p(:,:),gw(:),gamma,ak,theta
      real(dp),intent(in)::a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta
      integer::k,j,l,i,mjk,sticky,mtot,kocc,saux,ntot
      integer::ntr(size(p,1),size(p,2)),mbar(size(p,1),size(p,2))
      real(dp)::alpha,kappa,base(size(gw)),prow(size(gw)),eta,rate,prob,shape,nonsticky
      k=size(gw);call transition_counts(state,k,ntr);alpha=max((1.0_dp-theta)*ak,1.0e-8_dp);kappa=max(theta*ak,1.0e-8_dp)
      mbar=0;sticky=0;mtot=0
      do j=1,k;do l=1,k
         mjk=crt_draw(ntr(j,l),alpha*gw(l)+merge(kappa,0.0_dp,j==l));mtot=mtot+mjk
         if(j==l.and.mjk>0)then
            do i=1,mjk
               if(runif()<kappa/(kappa+alpha*gw(j)))then;sticky=sticky+1;else;mbar(j,l)=mbar(j,l)+1;end if
            end do
         else;mbar(j,l)=mjk;end if
      end do;end do
      do l=1,k;base(l)=gamma/real(k,dp)+real(sum(mbar(:,l)),dp);end do;call rdirichlet_rng(base,gw)
      alpha=max((1.0_dp-theta)*ak,1.0e-8_dp);kappa=max(theta*ak,1.0e-8_dp)
      do j=1,k
         do l=1,k;base(l)=alpha*gw(l)+real(ntr(j,l),dp)+merge(kappa,0.0_dp,j==l);end do
         call rdirichlet_rng(base,prow);p(j,:)=prow
      end do
      ! Concentration updates using Escobar-West/CRT auxiliaries.
      ntot=sum(mbar);kocc=count([(sum(mbar(:,l))>0,l=1,k)])
      if(ntot>0)then
         eta=rbeta(gamma+1.0_dp,real(ntot,dp));rate=b_gamma-log(max(eta,tiny(1.0_dp)))
         prob=(a_gamma+real(max(kocc,1)-1,dp))/(a_gamma+real(max(kocc,1)-1,dp)+real(ntot,dp)*rate)
         shape=a_gamma+real(max(kocc,1)-merge(0,1,runif()<prob),dp);gamma=rgamma_mt(max(shape,1.0e-6_dp),1.0_dp/rate)
      end if
      saux=0;rate=b_alpha
      do j=1,k
         ntot=sum(ntr(j,:));if(ntot<=0)cycle
         eta=rbeta(ak+1.0_dp,real(ntot,dp));rate=rate-log(max(eta,tiny(1.0_dp)))
         if(runif()<real(ntot,dp)/(real(ntot,dp)+ak))saux=saux+1
      end do
      shape=max(1.0e-6_dp,a_alpha+real(mtot-saux,dp));ak=rgamma_mt(shape,1.0_dp/max(rate,1.0e-8_dp))
      nonsticky=real(max(mtot-sticky,0),dp);theta=rbeta(a_theta+real(sticky,dp),b_theta+nonsticky)
   end subroutine sticky_hdp_update

   subroutine hdp_state_draw(y,x,beta,p,family,rho,state,prob,info)
      integer,intent(in)::y(:),family
      real(dp),intent(in)::x(:,:),beta(:,:),p(:,:),rho(:)
      integer,intent(inout)::state(:)
      real(dp),intent(out)::prob(:,:)
      integer,intent(out)::info
      integer::t,j,k
      real(dp)::le(size(y),size(beta,1)),pi0(size(beta,1)),mu
      k=size(beta,1);pi0=1.0_dp/real(k,dp)
      do t=1,size(y);do j=1,k
         mu=exp(max(-700.0_dp,min(700.0_dp,dot_product(x(t,:),beta(j,:)))))
         if(family==1)then;le(t,j)=poisson_logpmf(y(t),mu);else;le(t,j)=nb_logpmf(y(t),mu,rho(j));end if
      end do;end do
      call general_hmm_ffbs(le,p,pi0,state,prob,info)
   end subroutine hdp_state_draw

   subroutine initialize_hdp(p,gw)
      real(dp),intent(in)::p(:,:)
      real(dp),intent(out)::gw(size(p,1))
            integer::j
      gw=0.0_dp
      do j=1,size(p,1);gw=gw+p(j,:);end do
      gw=max(gw,1.0e-10_dp);gw=gw/sum(gw)
   end subroutine initialize_hdp

   function hdphmm_poisson(y,x,k,beta_start,p_start,gamma_start,ak_start,theta_start,b0,bprec, &
                           a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta,burnin,mcmc,thin,beta_step) result(res)
      integer,intent(in)::y(:),k,burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),p_start(:,:),gamma_start,ak_start,theta_start,b0(:),bprec(:,:)
      real(dp),intent(in)::a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta,beta_step
      type(hdp_count_result)::res
      integer::n,pdim,nstore,iter,keep,j,info
      real(dp)::beta(k,size(x,2)),p(k,k),gw(k),gamma,ak,theta,rhodummy(k)
      integer::state(size(y));real(dp)::prob(size(y),k)
      n=size(y);pdim=size(x,2);nstore=mcmc/thin
      if(size(x,1)/=n.or.any(y<0).or.any(shape(beta_start)/=[k,pdim]).or.any(shape(p_start)/=[k,k]).or. &
         size(b0)/=pdim.or.any(shape(bprec)/=[pdim,pdim]).or.gamma_start<=0.0_dp.or.ak_start<=0.0_dp.or. &
         theta_start<=0.0_dp.or.theta_start>=1.0_dp.or.thin<=0.or.nstore<=0)then;res%status=1;return;end if
      beta=beta_start;p=p_start;gamma=gamma_start;ak=ak_start;theta=theta_start;rhodummy=1.0_dp;state=1
      call initialize_hdp(p,gw)
      allocate(res%beta(nstore,k,pdim),res%p(nstore,k,k),res%state(nstore,n),res%state_prob(nstore,n,k), &
               res%gamma(nstore),res%alpha_kappa(nstore),res%theta(nstore));keep=0
      do iter=0,burnin+mcmc-1
         call hdp_state_draw(y,x,beta,p,1,rhodummy,state,prob,info);if(info/=0)then;res%status=10+info;return;end if
         do j=1,k;call beta_rw_update(y,x,state,j,beta(j,:),b0,bprec,1,1.0_dp,beta_step);end do
         call sticky_hdp_update(state,p,gw,gamma,ak,theta,a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta)
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;res%beta(keep,:,:)=beta;res%p(keep,:,:)=p;res%state(keep,:)=state;res%state_prob(keep,:,:)=prob
            res%gamma(keep)=gamma;res%alpha_kappa(keep)=ak;res%theta(keep)=theta
         end if
      end do
   end function hdphmm_poisson

   function hdphmm_negbin(y,x,k,beta_start,p_start,rho_start,gamma_start,ak_start,theta_start,b0,bprec, &
                          a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta,g,e,f,rho_step, &
                          burnin,mcmc,thin,beta_step) result(res)
      integer,intent(in)::y(:),k,burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),p_start(:,:),rho_start(:),gamma_start,ak_start,theta_start,b0(:),bprec(:,:)
      real(dp),intent(in)::a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta,g,e,f,rho_step(:),beta_step
      type(hdp_count_result)::res
      integer::n,pdim,nstore,iter,keep,j,info
      real(dp)::beta(k,size(x,2)),p(k,k),gw(k),rho(k),gamma,ak,theta
      integer::state(size(y));real(dp)::prob(size(y),k)
      n=size(y);pdim=size(x,2);nstore=mcmc/thin
      if(size(x,1)/=n.or.any(y<0).or.any(shape(beta_start)/=[k,pdim]).or.any(shape(p_start)/=[k,k]).or. &
         size(rho_start)/=k.or.size(rho_step)/=k.or.any(rho_start<=0.0_dp).or.size(b0)/=pdim.or. &
         any(shape(bprec)/=[pdim,pdim]).or.gamma_start<=0.0_dp.or.ak_start<=0.0_dp.or.theta_start<=0.0_dp.or. &
         theta_start>=1.0_dp.or.g<=0.0_dp.or.e<=0.0_dp.or.f<=0.0_dp.or.thin<=0.or.nstore<=0)then;res%status=1;return;end if
      beta=beta_start;p=p_start;rho=rho_start;gamma=gamma_start;ak=ak_start;theta=theta_start;state=1
      call initialize_hdp(p,gw)
      allocate(res%beta(nstore,k,pdim),res%p(nstore,k,k),res%state(nstore,n),res%state_prob(nstore,n,k),res%rho(nstore,k), &
               res%gamma(nstore),res%alpha_kappa(nstore),res%theta(nstore));keep=0
      do iter=0,burnin+mcmc-1
         call hdp_state_draw(y,x,beta,p,2,rho,state,prob,info);if(info/=0)then;res%status=10+info;return;end if
         do j=1,k
            call beta_rw_update(y,x,state,j,beta(j,:),b0,bprec,2,rho(j),beta_step)
            call rho_slice_update(y,x,state,j,beta(j,:),rho(j),rho_step(j),g,e,f)
         end do
         call sticky_hdp_update(state,p,gw,gamma,ak,theta,a_alpha,b_alpha,a_gamma,b_gamma,a_theta,b_theta)
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1
            res%beta(keep,:,:)=beta;res%p(keep,:,:)=p;res%state(keep,:)=state
            res%state_prob(keep,:,:)=prob;res%rho(keep,:)=rho
            res%gamma(keep)=gamma;res%alpha_kappa(keep)=ak;res%theta(keep)=theta
         end if
      end do
   end function hdphmm_negbin

   pure real(dp) function ztnb_dur_logpmf(d,r,omega) result(v)
      integer,intent(in)::d
      real(dp),intent(in)::r,omega
      real(dp)::z
      if(d<1.or.r<=0.0_dp.or.omega<=0.0_dp.or.omega>=1.0_dp)then;v=-huge(1.0_dp);return;end if
      z=1.0_dp-omega**r
      v=log_gamma(real(d,dp)+r)-log_gamma(r)-log_gamma(real(d+1,dp)) &
         +r*log(omega)+real(d,dp)*log(1.0_dp-omega)-log(max(z,tiny(1.0_dp)))
   end function ztnb_dur_logpmf

   subroutine hsmm_state_draw(y,x,beta,p,rho,omega,r,state,info)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:),beta(:,:),p(:,:),rho(:),omega(:),r
      integer,intent(out)::state(size(y));integer,intent(out)::info
      integer::n,k,t,j,d,u,l,prev,endt,startt,sel
      real(dp)::emit(size(y),size(beta,1)),cum(0:size(y),size(beta,1)),fwd(size(y),size(beta,1))
      real(dp)::pbar(size(p,1),size(p,2)),logs(size(p,1)),cand(size(y)*size(p,1)),mx,tot,uu,mu
      n=size(y);k=size(beta,1);info=0;pbar=p
      do j=1,k;pbar(j,j)=0.0_dp;if(sum(pbar(j,:))>0.0_dp)pbar(j,:)=pbar(j,:)/sum(pbar(j,:));end do
      do t=1,n;do j=1,k
         mu=exp(max(-700.0_dp,min(700.0_dp,dot_product(x(t,:),beta(j,:)))));emit(t,j)=nb_logpmf(y(t),mu,rho(j))
      end do;end do
      cum=0.0_dp;do j=1,k;do t=1,n;cum(t,j)=cum(t-1,j)+emit(t,j);end do;end do
      fwd=-huge(1.0_dp)
      do t=1,n;do j=1,k
         sel=0
         do d=1,t
            u=t-d+1
            if(u==1)then;logs(1)=-log(real(k,dp))+ztnb_dur_logpmf(d,r,omega(j))+cum(t,j)-cum(u-1,j);tot=logs(1)
            else
               do l=1,k
                  if(l==j.or.pbar(l,j)<=0.0_dp)then;logs(l)=-huge(1.0_dp);else;logs(l)=fwd(u-1,l)+log(pbar(l,j));end if
               end do
               tot=logsumexp(logs)+ztnb_dur_logpmf(d,r,omega(j))+cum(t,j)-cum(u-1,j)
            end if
            sel=sel+1;cand(sel)=tot
         end do
         fwd(t,j)=logsumexp(cand(1:sel))
      end do;end do
      state=1;endt=n
      ! Draw final state.
      logs=fwd(endt,:);mx=maxval(logs);tot=sum(exp(logs-mx));uu=runif()*tot;sel=1
      do j=1,k;uu=uu-exp(logs(j)-mx);if(uu<=0.0_dp)then;sel=j;exit;end if;end do; j=sel
      do while(endt>=1)
         sel=0
         do d=1,endt
            startt=endt-d+1
            if(startt==1)then;tot=-log(real(k,dp))+ztnb_dur_logpmf(d,r,omega(j))+cum(endt,j)-cum(startt-1,j)
            else
               do l=1,k
                  if(l==j.or.pbar(l,j)<=0.0_dp)then;logs(l)=-huge(1.0_dp);else;logs(l)=fwd(startt-1,l)+log(pbar(l,j));end if
               end do
               tot=logsumexp(logs)+ztnb_dur_logpmf(d,r,omega(j))+cum(endt,j)-cum(startt-1,j)
            end if
            sel=sel+1;cand(sel)=tot
         end do
         mx=maxval(cand(1:sel));tot=sum(exp(cand(1:sel)-mx));uu=runif()*tot;d=1
         do u=1,sel;uu=uu-exp(cand(u)-mx);if(uu<=0.0_dp)then;d=u;exit;end if;end do
         startt=endt-d+1;state(startt:endt)=j
         if(startt==1)exit
         do l=1,k
            if(l==j.or.pbar(l,j)<=0.0_dp)then;logs(l)=-huge(1.0_dp);else;logs(l)=fwd(startt-1,l)+log(pbar(l,j));end if
         end do
         mx=maxval(logs);tot=sum(exp(logs-mx));uu=runif()*tot;prev=1
         do l=1,k;uu=uu-exp(logs(l)-mx);if(uu<=0.0_dp)then;prev=l;exit;end if;end do
         j=prev;endt=startt-1
      end do
   end subroutine hsmm_state_draw

   function hdphsmm_negbin(y,x,k,beta_start,p_start,rho_start,omega_start,gamma_start,alpha_start,b0,bprec, &
                           a_alpha,b_alpha,a_gamma,b_gamma,a_omega,b_omega,g,e,f,dur_r,rho_step, &
                           burnin,mcmc,thin,beta_step) result(res)
      integer,intent(in)::y(:),k,burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),p_start(:,:),rho_start(:),omega_start(:),gamma_start,alpha_start,b0(:),bprec(:,:)
      real(dp),intent(in)::a_alpha,b_alpha,a_gamma,b_gamma,a_omega,b_omega,g,e,f,dur_r,rho_step(:),beta_step
      type(hdp_count_result)::res
      integer::n,pdim,nstore,iter,keep,j,info,t,nseg,totdur
      real(dp)::beta(k,size(x,2)),p(k,k),gw(k),rho(k),omega(k),gamma,ak,theta
      integer::state(size(y));real(dp)::prob(size(y),k)
      n=size(y);pdim=size(x,2);nstore=mcmc/thin
      if(size(x,1)/=n.or.any(y<0).or.any(shape(beta_start)/=[k,pdim]).or.any(shape(p_start)/=[k,k]).or. &
         size(rho_start)/=k.or.size(omega_start)/=k.or.size(rho_step)/=k.or.any(rho_start<=0.0_dp).or.any(omega_start<=0.0_dp).or. &
         any(omega_start>=1.0_dp).or.gamma_start<=0.0_dp.or.alpha_start<=0.0_dp.or.a_omega<=0.0_dp.or.b_omega<=0.0_dp.or. &
         dur_r<=0.0_dp.or.thin<=0.or.nstore<=0)then;res%status=1;return;end if
      beta=beta_start;p=p_start;rho=rho_start;omega=omega_start;gamma=gamma_start;ak=alpha_start;theta=1.0e-6_dp;state=1
      call initialize_hdp(p,gw)
      allocate(res%beta(nstore,k,pdim),res%p(nstore,k,k),res%state(nstore,n), &
               res%state_prob(nstore,n,k),res%rho(nstore,k),res%omega(nstore,k), &
               res%gamma(nstore),res%alpha_kappa(nstore),res%theta(nstore));keep=0
      do iter=0,burnin+mcmc-1
         call hsmm_state_draw(y,x,beta,p,rho,omega,dur_r,state,info);if(info/=0)then;res%status=10+info;return;end if
         do j=1,k
            call beta_rw_update(y,x,state,j,beta(j,:),b0,bprec,2,rho(j),beta_step)
            call rho_slice_update(y,x,state,j,beta(j,:),rho(j),rho_step(j),g,e,f)
            nseg=0;totdur=0;t=1
            do while(t<=n)
               if(state(t)==j)then
                  nseg=nseg+1
                  do
                     if(t>n)exit
                     if(state(t)/=j)exit
                     totdur=totdur+1;t=t+1
                  end do
               else
                  t=t+1
               end if
            end do
            omega(j)=rbeta(a_omega+dur_r*real(nseg,dp),b_omega+real(totdur,dp))
         end do
         ! HSMM has no sticky self transition. Draw rows from alpha*global weights + segment transitions.
         call hsmm_transition_update(state,p,gw,gamma,ak,a_alpha,b_alpha,a_gamma,b_gamma)
         prob=0.0_dp;do t=1,n;prob(t,state(t))=1.0_dp;end do
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1
            res%beta(keep,:,:)=beta;res%p(keep,:,:)=p;res%state(keep,:)=state
            res%state_prob(keep,:,:)=prob;res%rho(keep,:)=rho;res%omega(keep,:)=omega
            res%gamma(keep)=gamma;res%alpha_kappa(keep)=ak;res%theta(keep)=0.0_dp
         end if
      end do
   end function hdphsmm_negbin

   subroutine hsmm_transition_update(state,p,gw,gamma,alpha,a_alpha,b_alpha,a_gamma,b_gamma)
      integer,intent(in)::state(:)
      real(dp),intent(inout)::p(:,:),gw(:),gamma,alpha
      real(dp),intent(in)::a_alpha,b_alpha,a_gamma,b_gamma
      integer::k,j,l,t,mjk,mtot,kocc,ntot,saux
      integer::ntr(size(p,1),size(p,2)),m(size(p,1),size(p,2))
      real(dp)::base(size(gw)),prow(size(gw)),eta,rate,prob,shape
      k=size(gw);ntr=0
      do t=1,size(state)-1
         if(state(t)/=state(t+1))ntr(state(t),state(t+1))=ntr(state(t),state(t+1))+1
      end do
      m=0;mtot=0
      do j=1,k;do l=1,k
         if(j/=l)then;mjk=crt_draw(ntr(j,l),max(alpha*gw(l),1.0e-10_dp));m(j,l)=mjk;mtot=mtot+mjk;end if
      end do;end do
      do l=1,k;base(l)=gamma/real(k,dp)+real(sum(m(:,l)),dp);end do;call rdirichlet_rng(base,gw)
      do j=1,k
         do l=1,k;base(l)=merge(1.0e-12_dp,alpha*gw(l)+real(ntr(j,l),dp),j==l);end do
         call rdirichlet_rng(base,prow);p(j,:)=prow;p(j,j)=0.0_dp;p(j,:)=p(j,:)/max(sum(p(j,:)),tiny(1.0_dp))
      end do
      kocc=count([(sum(m(:,l))>0,l=1,k)]);if(mtot>0)then
         eta=rbeta(gamma+1.0_dp,real(mtot,dp));rate=b_gamma-log(max(eta,tiny(1.0_dp)))
         prob=(a_gamma+real(max(kocc,1)-1,dp))/(a_gamma+real(max(kocc,1)-1,dp)+real(mtot,dp)*rate)
         shape=a_gamma+real(max(kocc,1)-merge(0,1,runif()<prob),dp);gamma=rgamma_mt(max(shape,1.0e-6_dp),1.0_dp/rate)
      end if
      saux=0;rate=b_alpha
      do j=1,k;ntot=sum(ntr(j,:));if(ntot>0)then
         eta=rbeta(alpha+1.0_dp,real(ntot,dp))
         rate=rate-log(max(eta,tiny(1.0_dp)))
         if(runif()<real(ntot,dp)/(real(ntot,dp)+alpha))saux=saux+1
      end if;end do
      alpha=rgamma_mt(max(a_alpha+real(mtot-saux,dp),1.0e-6_dp),1.0_dp/max(rate,1.0e-8_dp))
   end subroutine hsmm_transition_update
end module mcmcpack_hdp_hmm
