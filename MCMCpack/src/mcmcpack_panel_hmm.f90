! SPDX-License-Identifier: GPL-3.0-only
! Gaussian panel hidden-Markov samplers translated from cHMMpanelFE.cc/cHMMpanelRE.cc.
! Data are supplied as dense subject x time arrays; this replaces R's long-form indexing glue.
module mcmcpack_panel_hmm
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : rnorm, rbeta, rmvnorm, rinvgamma_rng
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_distributions, only : riwish
   use mcmcpack_special_utils, only : normal_logpdf, general_hmm_ffbs
   implicit none
   private
   public :: panel_fe_result, panel_re_result, hmm_panel_fe, hmm_panel_re

   type :: panel_fe_result
      real(dp),allocatable :: beta(:,:)
      real(dp),allocatable :: delta(:,:,:)
      real(dp),allocatable :: sigma2(:,:,:)
      integer,allocatable :: state(:,:,:)
      real(dp),allocatable :: state_prob(:,:,:,:)
      integer :: status=0
   end type panel_fe_result
   type :: panel_re_result
      real(dp),allocatable :: beta(:,:,:)
      real(dp),allocatable :: sigma2(:,:)
      real(dp),allocatable :: d(:,:,:,:)
      integer,allocatable :: state(:,:)
      real(dp),allocatable :: state_prob(:,:,:)
      integer :: status=0
   end type panel_re_result
contains
   subroutine ordered_transition_draw(state,s,p0,p)
      integer,intent(in)::state(:),s
      real(dp),intent(in)::p0(:,:)
      real(dp),intent(out)::p(s,s)
      integer::j,t,nsame,nnext
      real(dp)::v
      p=0.0_dp
      if(s==1)then;p(1,1)=1.0_dp;return;end if
      do j=1,s-1
         nsame=0;nnext=0
         do t=1,size(state)-1
            if(state(t)==j.and.state(t+1)==j)nsame=nsame+1
            if(state(t)==j.and.state(t+1)==j+1)nnext=nnext+1
         end do
         v=rbeta(max(p0(j,j),1.0e-6_dp)+real(nsame,dp),max(p0(j,j+1),1.0e-6_dp)+real(nnext,dp))
         p(j,j)=v;p(j,j+1)=1.0_dp-v
      end do
      p(s,s)=1.0_dp
   end subroutine ordered_transition_draw

   function hmm_panel_fe(y,x,nstate,beta_start,delta_start,sigma2_start,b0,b0prec, &
                         delta0,delta_prec,c0,d0,p_start,p0,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:,:),x(:,:,:),beta_start(:),delta_start(:,:),sigma2_start(:,:)
      integer,intent(in)::nstate(:),burnin,mcmc,thin
      real(dp),intent(in)::b0(:),b0prec(:,:),delta0,delta_prec,c0,d0,p_start(:,:,:),p0(:,:,:)
      type(panel_fe_result)::res
      integer::ns,nt,k,smax,nstore,iter,keep,i,t,j,st,info
      real(dp)::beta(size(beta_start)),delta(size(delta_start,1),size(delta_start,2))
      real(dp)::sig2(size(sigma2_start,1),size(sigma2_start,2)),p(size(p_start,1),size(p_start,2),size(p_start,3))
      integer,allocatable::state(:,:),path(:)
      real(dp),allocatable::prob(:,:,:),le(:,:),pt(:,:),pi0(:),pprob(:,:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:)
      real(dp)::v,m,sse,r
      ns=size(y,1);nt=size(y,2);k=size(x,3);smax=maxval(nstate);nstore=mcmc/thin
      if(size(x,1)/=ns.or.size(x,2)/=nt.or.size(beta_start)/=k.or.size(nstate)/=ns.or. &
         any(nstate<1).or.any(nstate>smax).or.any(shape(delta_start)/=[ns,smax]).or. &
         any(shape(sigma2_start)/=[ns,smax]).or.any(shape(p_start)/=[ns,smax,smax]).or. &
         any(shape(p0)/=[ns,smax,smax]).or.size(b0)/=k.or.any(shape(b0prec)/=[k,k]).or. &
         delta_prec<0.0_dp.or.c0<=0.0_dp.or.d0<=0.0_dp.or.thin<=0.or.nstore<=0)then;res%status=1;return;end if
      beta=beta_start;delta=delta_start;sig2=sigma2_start;p=p_start
      allocate(state(ns,nt),prob(ns,nt,smax));state=1;prob=0.0_dp
      allocate(res%beta(nstore,k),res%delta(nstore,ns,smax),res%sigma2(nstore,ns,smax), &
               res%state(nstore,ns,nt),res%state_prob(nstore,ns,nt,smax))
      allocate(prec(k,k),cov(k,k),rhs(k),mu(k),draw(k));keep=0
      do iter=0,burnin+mcmc-1
         ! Subject-specific ordered state paths.
         do i=1,ns
            allocate(le(nt,nstate(i)),pt(nstate(i),nstate(i)),pi0(nstate(i)),path(nt),pprob(nt,nstate(i)))
            pt=p(i,1:nstate(i),1:nstate(i));pi0=0.0_dp;pi0(1)=1.0_dp
            do t=1,nt;do j=1,nstate(i)
               r=y(i,t)-dot_product(x(i,t,:),beta)-delta(i,j)
               le(t,j)=normal_logpdf(r,0.0_dp,sqrt(max(sig2(i,j),1.0e-12_dp)))
            end do;end do
            call general_hmm_ffbs(le,pt,pi0,path,pprob,info);state(i,:)=path;prob(i,:,1:nstate(i))=pprob
            if(info/=0)then;res%status=10+info;return;end if
            call ordered_transition_draw(state(i,:),nstate(i),p0(i,1:nstate(i),1:nstate(i)),pt)
            p(i,1:nstate(i),1:nstate(i))=pt
            deallocate(le,pt,pi0,path,pprob)
         end do
         ! State intercepts and variances.
         do i=1,ns;do j=1,nstate(i)
            v=delta_prec;m=delta_prec*delta0
            do t=1,nt
               if(state(i,t)/=j)cycle
               v=v+1.0_dp/sig2(i,j)
               m=m+(y(i,t)-dot_product(x(i,t,:),beta))/sig2(i,j)
            end do
            if(v>0.0_dp)delta(i,j)=rnorm(m/v,sqrt(1.0_dp/v))
            sse=0.0_dp;st=0
            do t=1,nt
               if(state(i,t)==j)then;r=y(i,t)-dot_product(x(i,t,:),beta)-delta(i,j);sse=sse+r*r;st=st+1;end if
            end do
            sig2(i,j)=rinvgamma_rng(0.5_dp*(c0+real(st,dp)),0.5_dp*(d0+sse))
         end do;end do
         ! Common fixed-effect beta.
         prec=b0prec;rhs=matmul(b0prec,b0)
         do i=1,ns;do t=1,nt
            j=state(i,t);r=1.0_dp/max(sig2(i,j),1.0e-12_dp)
            prec=prec+r*spread(x(i,t,:),2,k)*spread(x(i,t,:),1,k)
            rhs=rhs+r*x(i,t,:)*(y(i,t)-delta(i,j))
         end do;end do
         call inv_spd(prec,cov,info);if(info/=0)then;res%status=20+info;return;end if
         mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=30+info;return;end if;beta=draw
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;res%beta(keep,:)=beta;res%delta(keep,:,:)=delta;res%sigma2(keep,:,:)=sig2
            res%state(keep,:,:)=state;res%state_prob(keep,:,:,:)=prob
         end if
      end do
   end function hmm_panel_fe

   function hmm_panel_re(y,x,w,nstate,beta_start,sigma2_start,d_start,b0,b0prec,c0,d0,r0,rscale, &
                         p_start,p0,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:,:),x(:,:,:),w(:,:,:),beta_start(:,:),sigma2_start(:),d_start(:,:,:)
      integer,intent(in)::nstate,burnin,mcmc,thin
      real(dp),intent(in)::b0(:),b0prec(:,:),c0,d0,r0,rscale(:,:),p_start(:,:),p0(:,:)
      type(panel_re_result)::res
      integer::ns,nt,k,q,s,nstore,iter,keep,i,t,j,a,b,info,nj
      real(dp)::beta(size(beta_start,1),size(beta_start,2)),sig2(size(sigma2_start))
      real(dp)::d(size(d_start,1),size(d_start,2),size(d_start,3))
      real(dp)::p(size(p_start,1),size(p_start,2)),le(size(y,2),nstate),pi0(nstate),prob(size(y,2),nstate)
      integer::state(size(y,2))
      real(dp),allocatable::bi(:,:,:),dinv(:,:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:),scale(:,:)
      real(dp)::rr,sse
      ns=size(y,1);nt=size(y,2);k=size(x,3);q=size(w,3);s=nstate;nstore=mcmc/thin
      if(size(x,1)/=ns.or.size(x,2)/=nt.or.size(w,1)/=ns.or.size(w,2)/=nt.or. &
         any(shape(beta_start)/=[k,s]).or.size(sigma2_start)/=s.or.any(shape(d_start)/=[q,q,s]).or. &
         size(b0)/=k.or.any(shape(b0prec)/=[k,k]).or.any(shape(rscale)/=[q,q]).or. &
         any(shape(p_start)/=[s,s]).or.any(shape(p0)/=[s,s]).or.c0<=0.0_dp.or.d0<=0.0_dp.or. &
         r0<real(q,dp).or.thin<=0.or.nstore<=0)then;res%status=1;return;end if
      beta=beta_start;sig2=sigma2_start;d=d_start;p=p_start;state=1;pi0=0.0_dp;pi0(1)=1.0_dp
      allocate(bi(ns,q,s));bi=0.0_dp
      allocate(res%beta(nstore,k,s),res%sigma2(nstore,s),res%d(nstore,q,q,s), &
               res%state(nstore,nt),res%state_prob(nstore,nt,s))
      keep=0
      do iter=0,burnin+mcmc-1
         ! Common time-regime path conditional on current random effects.
         do t=1,nt;do j=1,s
            le(t,j)=0.0_dp
            do i=1,ns
               rr=y(i,t)-dot_product(x(i,t,:),beta(:,j))-dot_product(w(i,t,:),bi(i,:,j))
               le(t,j)=le(t,j)+normal_logpdf(rr,0.0_dp,sqrt(max(sig2(j),1.0e-12_dp)))
            end do
         end do;end do
         call general_hmm_ffbs(le,p,pi0,state,prob,info);if(info/=0)then;res%status=10+info;return;end if
         call ordered_transition_draw(state,s,p0,p)
         do j=1,s
            allocate(dinv(q,q))
            call inv_spd(d(:,:,j),dinv,info);if(info/=0)then;res%status=20+info;return;end if
            ! Random effects by subject and regime.
            do i=1,ns
               allocate(prec(q,q),cov(q,q),rhs(q),mu(q),draw(q));prec=dinv;rhs=0.0_dp
               do t=1,nt
                  if(state(t)/=j)cycle
                  rr=1.0_dp/max(sig2(j),1.0e-12_dp)
                  prec=prec+rr*spread(w(i,t,:),2,q)*spread(w(i,t,:),1,q)
                  rhs=rhs+rr*w(i,t,:)*(y(i,t)-dot_product(x(i,t,:),beta(:,j)))
               end do
               call inv_spd(prec,cov,info);if(info/=0)then;res%status=30+info;return;end if
               mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=40+info;return;end if;bi(i,:,j)=draw
               deallocate(prec,cov,rhs,mu,draw)
            end do
            ! Fixed effects for this regime.
            allocate(prec(k,k),cov(k,k),rhs(k),mu(k),draw(k));prec=b0prec;rhs=matmul(b0prec,b0)
            do i=1,ns;do t=1,nt
               if(state(t)/=j)cycle;rr=1.0_dp/max(sig2(j),1.0e-12_dp)
               prec=prec+rr*spread(x(i,t,:),2,k)*spread(x(i,t,:),1,k)
               rhs=rhs+rr*x(i,t,:)*(y(i,t)-dot_product(w(i,t,:),bi(i,:,j)))
            end do;end do
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=50+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=60+info;return;end if;beta(:,j)=draw
            deallocate(prec,cov,rhs,mu,draw)
            ! Residual variance.
            sse=0.0_dp;nj=0
            do i=1,ns;do t=1,nt
               if(state(t)==j)then
                  rr=y(i,t)-dot_product(x(i,t,:),beta(:,j)) &
                     -dot_product(w(i,t,:),bi(i,:,j))
                  sse=sse+rr*rr;nj=nj+1
               end if
            end do;end do
            sig2(j)=rinvgamma_rng(0.5_dp*(c0+real(nj,dp)),0.5_dp*(d0+sse))
            ! Random-effect covariance D_j ~ IW(r0+n, R0+sum b_i b_i').
            allocate(scale(q,q));scale=rscale
            do i=1,ns
               do a=1,q;do b=1,q;scale(a,b)=scale(a,b)+bi(i,a,j)*bi(i,b,j);end do;end do
            end do
            call riwish(r0+real(ns,dp),scale,d(:,:,j),info);deallocate(scale,dinv)
            if(info/=0)then;res%status=70+info;return;end if
         end do
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;res%beta(keep,:,:)=beta;res%sigma2(keep,:)=sig2;res%d(keep,:,:,:)=d
            res%state(keep,:)=state;res%state_prob(keep,:,:)=prob
         end if
      end do
   end function hmm_panel_re
end module mcmcpack_panel_hmm
