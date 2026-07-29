! SPDX-License-Identifier: GPL-3.0-only
module nmof_portfolio
   use nmof_kinds, only: dp
   use nmof_linalg, only: solve_linear, covariance_matrix, column_means
   use nmof_qp, only: solve_qp_active_set, project_feasible, project_budget_box
   use nmof_types, only: frontier_result, nmof_ok, nmof_invalid_input, nmof_numerical_failure
   implicit none
   private
   public :: equal_risk_contribution, minimum_variance, mean_variance_portfolio
   public :: mean_variance_frontier, maximum_sharpe, tracking_portfolio
   public :: minimum_cvar, minimum_mad
contains
   subroutine minimum_variance(covariance,w,lower,upper,groups,group_lower,group_upper,status)
      real(dp),intent(in)::covariance(:,:)
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::lower(:),upper(:),groups(:,:),group_lower(:),group_upper(:)
      integer,intent(out),optional::status
      real(dp),allocatable::lo(:),hi(:),c(:,:),cv(:),d(:,:),dv(:),q(:,:),qv(:),x0(:),feasible(:)
      integer::n,info,it
      n=size(w)
      if(any(shape(covariance)/=[n,n]).or.n<1) then; if(present(status)) status=nmof_invalid_input; return; end if
      call default_bounds(n,lower,upper,lo,hi)
      call build_constraints(n,lo,hi,groups,group_lower,group_upper,.true.,c,cv,d,dv)
      allocate(q(n,n),qv(n),x0(n)); q=2.0_dp*covariance; call regularize(q); qv=0.0_dp; x0=1.0_dp/real(n,dp)
      allocate(feasible(size(x0))); call project_feasible(x0,c,cv,d,dv,feasible,info); x0=feasible
      if(info/=nmof_ok) then; if(present(status)) status=info; return; end if
      call solve_qp_active_set(q,qv,c,cv,d,dv,x0,w,info,it)
      if(present(status)) status=info
   end subroutine minimum_variance

   subroutine mean_variance_portfolio(mean_returns,covariance,min_return,w,lower,upper,lambda,groups,group_lower,group_upper,status)
      real(dp),intent(in)::mean_returns(:),covariance(:,:),min_return
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::lower(:),upper(:),lambda(:),groups(:,:),group_lower(:),group_upper(:)
      integer,intent(out),optional::status
      real(dp),allocatable::lo(:),hi(:),c(:,:),cv(:),d(:,:),dv(:),d2(:,:),dv2(:),q(:,:),qv(:),x0(:),feasible(:)
      real(dp)::l1,l2
      integer::n,info,it,m
      n=size(w)
      if(size(mean_returns)/=n.or.any(shape(covariance)/=[n,n])) then; if(present(status)) status=nmof_invalid_input; return; end if
      call default_bounds(n,lower,upper,lo,hi); call build_constraints(n,lo,hi,groups,group_lower,group_upper,.true.,c,cv,d,dv)
      l1=0.0_dp; l2=1.0_dp
      if(present(lambda)) then
         if(size(lambda)==1) then; l1=lambda(1); l2=1.0_dp-l1; else; l1=lambda(1); l2=lambda(2); end if
      else
         m=size(d,1); allocate(d2(m+1,n),dv2(m+1)); if(m>0) then; d2(1:m,:)=d; dv2(1:m)=dv; end if
         d2(m+1,:)=-mean_returns; dv2(m+1)=-min_return; call move_alloc(d2,d); call move_alloc(dv2,dv)
      end if
      allocate(q(n,n),qv(n),x0(n)); q=2.0_dp*l2*covariance; call regularize(q); qv=-l1*mean_returns; x0=1.0_dp/real(n,dp)
      allocate(feasible(size(x0))); call project_feasible(x0,c,cv,d,dv,feasible,info); x0=feasible
      if(info/=nmof_ok) then; if(present(status)) status=info; return; end if
      call solve_qp_active_set(q,qv,c,cv,d,dv,x0,w,info,it)
      if(present(status)) status=info
   end subroutine mean_variance_portfolio

   function mean_variance_frontier(mean_returns,covariance,n_points,lower,upper,risk_free,groups,group_lower,group_upper) result(frontier)
      real(dp),intent(in)::mean_returns(:),covariance(:,:)
      integer,intent(in),optional::n_points
      real(dp),intent(in),optional::lower(:),upper(:),risk_free,groups(:,:),group_lower(:),group_upper(:)
      type(frontier_result)::frontier
      real(dp),allocatable::w(:),lam(:),lo(:),hi(:)
      real(dp)::rf,target
      integer::n,np,i,st
      n=size(mean_returns); np=50; if(present(n_points)) np=n_points
      if(any(shape(covariance)/=[n,n]).or.np<1) then; frontier%status=nmof_invalid_input; return; end if
      allocate(frontier%returns(np),frontier%volatility(np),frontier%portfolios(n+merge(1,0,present(risk_free)),np),w(n),lam(1))
      call default_bounds(n,lower,upper,lo,hi)
      if(.not.present(risk_free)) then
         do i=1,np
            lam(1)=0.0001_dp+0.9998_dp*real(i-1,dp)/real(max(1,np-1),dp)
            call mean_variance_portfolio(mean_returns,covariance,0.0_dp,w,lo,hi,lam,groups,group_lower,group_upper,st)
            if(st/=nmof_ok) then; frontier%status=st; return; end if
            frontier%returns(i)=dot_product(mean_returns,w); frontier%volatility(i)=sqrt(max(0.0_dp,dot_product(w,matmul(covariance,w))))
            frontier%portfolios(:,i)=w
         end do
      else
         rf=risk_free
         do i=1,np
            target=rf+(maxval(mean_returns)-rf)*real(i-1,dp)/real(max(1,np-1),dp)
            call mean_variance_portfolio(mean_returns,covariance,target,w,lo,hi,groups=groups,group_lower=group_lower,group_upper=group_upper,status=st)
            if(st/=nmof_ok) then; frontier%status=st; return; end if
            frontier%portfolios(1:n,i)=w; frontier%portfolios(n+1,i)=1.0_dp-sum(w)
            frontier%returns(i)=dot_product(mean_returns,w)+(1.0_dp-sum(w))*rf
            frontier%volatility(i)=sqrt(max(0.0_dp,dot_product(w,matmul(covariance,w))))
         end do
      end if
      frontier%status=nmof_ok
   end function mean_variance_frontier

   subroutine maximum_sharpe(mean_returns,covariance,w,min_return,lower,upper,status)
      real(dp),intent(in)::mean_returns(:),covariance(:,:)
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::min_return,lower(:),upper(:)
      integer,intent(out),optional::status
      real(dp),allocatable::y(:),c(:,:),cv(:),d(:,:),dv(:),q(:,:),qv(:),x0(:),feasible(:),lo(:),hi(:)
      real(dp)::mr
      integer::n,i,info,it,m
      n=size(w); mr=1.0_dp; if(present(min_return)) mr=min_return
      if(any(shape(covariance)/=[n,n]).or.size(mean_returns)/=n) then; if(present(status)) status=nmof_invalid_input; return; end if
      if(.not.present(lower).and..not.present(upper)) then
         allocate(y(n)); call solve_linear(covariance,mean_returns,y,info)
         if(info==0.and.abs(sum(y))>tiny(1.0_dp)) then; w=y/sum(y); if(present(status)) status=nmof_ok; else; if(present(status)) status=nmof_numerical_failure; end if
         return
      end if
      call default_unbounded(n,lower,upper,lo,hi)
      allocate(c(1,n),cv(1)); c(1,:)=mean_returns; cv(1)=mr
      m=count(hi<huge(1.0_dp)/4)+count(lo>-huge(1.0_dp)/4); allocate(d(m,n),dv(m)); m=0
      do i=1,n
         if(hi(i)<huge(1.0_dp)/4) then; m=m+1; d(m,:)=-hi(i); d(m,i)=d(m,i)+1.0_dp; dv(m)=0.0_dp; end if
         if(lo(i)>-huge(1.0_dp)/4) then; m=m+1; d(m,:)=lo(i); d(m,i)=d(m,i)-1.0_dp; dv(m)=0.0_dp; end if
      end do
      allocate(q(n,n),qv(n),x0(n),y(n)); q=covariance; call regularize(q); qv=0.0_dp; x0=mr*mean_returns/max(dot_product(mean_returns,mean_returns),tiny(1.0_dp))
      allocate(feasible(size(x0))); call project_feasible(x0,c,cv,d,dv,feasible,info); x0=feasible
      if(info/=nmof_ok) then; if(present(status)) status=info; return; end if
      call solve_qp_active_set(q,qv,c,cv,d,dv,x0,y,info,it)
      if(abs(sum(y))>tiny(1.0_dp)) w=y/sum(y)
      if(present(status)) status=info
   end subroutine maximum_sharpe

   subroutine tracking_portfolio(covariance,w,lower,upper,returns,objective,status)
      real(dp),intent(in),optional::covariance(:,:),returns(:,:)
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::lower(:),upper(:)
      character(len=*),intent(in),optional::objective
      integer,intent(out),optional::status
      real(dp),allocatable::cov(:,:),q(:,:),qv(:),lo(:),hi(:),c(:,:),cv(:),d(:,:),dv(:),x0(:),feasible(:)
      integer::n,info,it
      character(len=24)::obj
      n=size(w); obj='variance'; if(present(objective)) obj=lowercase(objective)
      call default_bounds(n,lower,upper,lo,hi); call build_constraints(n,lo,hi,budget=.true.,c=c,cv=cv,d=d,dv=dv)
      allocate(q(n,n),qv(n),x0(n))
      if(index(obj,'sum')==1) then
         if(.not.present(returns).or.size(returns,2)/=n+1) then; if(present(status)) status=nmof_invalid_input; return; end if
         q=matmul(transpose(returns(:,2:)),returns(:,2:)); qv=-matmul(transpose(returns(:,2:)),returns(:,1))
      else
         if(present(covariance)) then; cov=covariance
         else if(present(returns)) then; cov=covariance_matrix(returns)
         else; if(present(status)) status=nmof_invalid_input; return; end if
         if(any(shape(cov)/=[n+1,n+1])) then; if(present(status)) status=nmof_invalid_input; return; end if
         q=cov(2:,2:); qv=-cov(1,2:)
      end if
      q=2.0_dp*q; qv=2.0_dp*qv; call regularize(q); x0=1.0_dp/real(n,dp); allocate(feasible(size(x0))); call project_feasible(x0,c,cv,d,dv,feasible,info); x0=feasible
      if(info==nmof_ok) call solve_qp_active_set(q,qv,c,cv,d,dv,x0,w,info,it)
      if(present(status)) status=info
   end subroutine tracking_portfolio

   subroutine equal_risk_contribution(covariance,w,lower,upper,status,tol,maxiter)
      real(dp),intent(in)::covariance(:,:)
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::lower(:),upper(:),tol
      integer,intent(out),optional::status
      integer,intent(in),optional::maxiter
      real(dp),allocatable::x(:),sx(:),lo(:),hi(:),rc(:)
      real(dp)::eps,aux,xnew,diff,scale
      integer::n,i,k,mit,info
      n=size(w); eps=1e-8_dp; if(present(tol)) eps=tol; mit=1000; if(present(maxiter)) mit=maxiter
      if(any(shape(covariance)/=[n,n]).or.any([(covariance(i,i),i=1,n)]<=0.0_dp)) then; if(present(status)) status=nmof_invalid_input; return; end if
      allocate(x(n),sx(n),rc(n)); scale=sum(covariance); if(scale<=0) scale=sum([(covariance(i,i),i=1,n)])
      x=sqrt(1.0_dp/scale); sx=matmul(covariance,x)
      do k=1,mit
         do i=1,n
            aux=x(i)*covariance(i,i)-sx(i); xnew=(aux+sqrt(aux*aux+4.0_dp*covariance(i,i)/real(n,dp)))/(2.0_dp*covariance(i,i))
            diff=xnew-x(i); sx=sx+covariance(:,i)*diff; x(i)=xnew
         end do
         w=x/sum(x); rc=w*matmul(covariance,w); rc=rc/sum(rc)
         if(maxval(abs(rc-1.0_dp/real(n,dp)))<eps) exit
      end do
      if(present(lower).or.present(upper)) then
         call default_bounds(n,lower,upper,lo,hi); x=w; call project_budget_box(x,lo,hi,w,info)
      end if
      if(present(status)) status=merge(nmof_ok,nmof_numerical_failure,k<=mit)
   end subroutine equal_risk_contribution

   subroutine minimum_cvar(returns,w,q,lower,upper,min_return,mean_returns,groups,group_lower,group_upper,status,maxiter,tol)
      real(dp),intent(in)::returns(:,:)
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::q,lower(:),upper(:),min_return,mean_returns(:),groups(:,:),group_lower(:),group_upper(:),tol
      integer,intent(out),optional::status
      integer,intent(in),optional::maxiter
      real(dp),allocatable::lo(:),hi(:),c(:,:),cv(:),d(:,:),dv(:),d2(:,:),dv2(:),x0(:),grad(:),loss(:),mu(:)
      real(dp)::tail,alpha,ga,step,eps,old,newobj
      integer::n,ns,it,mit,info,m,count_tail
      n=size(w); ns=size(returns,1); tail=0.1_dp; if(present(q)) tail=q; mit=20000; if(present(maxiter)) mit=maxiter; eps=1e-8_dp; if(present(tol)) eps=tol
      if(size(returns,2)/=n.or.tail<=0.or.tail>=1) then; if(present(status)) status=nmof_invalid_input; return; end if
      call default_bounds(n,lower,upper,lo,hi); call build_constraints(n,lo,hi,groups,group_lower,group_upper,.true.,c,cv,d,dv)
      if(present(min_return)) then
         if(present(mean_returns)) then; mu=mean_returns; else; mu=column_means(returns); end if
         m=size(d,1); allocate(d2(m+1,n),dv2(m+1)); if(m>0) then; d2(1:m,:)=d; dv2(1:m)=dv; end if
         d2(m+1,:)=-mu; dv2(m+1)=-min_return; call move_alloc(d2,d); call move_alloc(dv2,dv)
      end if
      allocate(x0(n),grad(n),loss(ns)); x0=1.0_dp/real(n,dp); call project_feasible(x0,c,cv,d,dv,w,info)
      if(info/=nmof_ok) then; if(present(status)) status=info; return; end if
      loss=-matmul(returns,w); alpha=quantile_simple(loss,1.0_dp-tail); old=huge(1.0_dp)
      do it=1,mit
         loss=-matmul(returns,w); count_tail=count(loss>alpha)
         grad=0.0_dp; if(count_tail>0) grad=-sum(returns,dim=1,mask=spread(loss>alpha,2,n))/(tail*real(ns,dp))
         ga=1.0_dp-real(count_tail,dp)/(tail*real(ns,dp)); step=0.05_dp/sqrt(real(it,dp))
         x0=w-step*grad; call project_feasible(x0,c,cv,d,dv,w,info); alpha=alpha-step*ga
         newobj=alpha+sum(max(loss-alpha,0.0_dp))/(tail*real(ns,dp))
         if(abs(newobj-old)<eps) exit; old=newobj
      end do
      if(present(status)) status=merge(nmof_ok,nmof_numerical_failure,it<=mit)
   contains
      function quantile_simple(x,p) result(v)
         use nmof_math, only: quantile_type7
         real(dp),intent(in)::x(:),p; real(dp)::v; v=quantile_type7(x,p)
      end function quantile_simple
   end subroutine minimum_cvar

   subroutine minimum_mad(returns,w,lower,upper,min_return,mean_returns,demean,status,maxiter,tol)
      real(dp),intent(in)::returns(:,:)
      real(dp),intent(out)::w(:)
      real(dp),intent(in),optional::lower(:),upper(:),min_return,mean_returns(:),tol
      logical,intent(in),optional::demean
      integer,intent(out),optional::status
      integer,intent(in),optional::maxiter
      real(dp),allocatable::lo(:),hi(:),c(:,:),cv(:),d(:,:),dv(:),d2(:,:),dv2(:),x0(:),grad(:),r0(:,:),rw(:),mu(:)
      real(dp)::step,old,obj,eps
      integer::n,ns,it,mit,info,m
      logical::dm
      n=size(w); ns=size(returns,1); dm=.true.; if(present(demean)) dm=demean; mit=20000; if(present(maxiter)) mit=maxiter; eps=1e-8_dp; if(present(tol)) eps=tol
      if(size(returns,2)/=n) then; if(present(status)) status=nmof_invalid_input; return; end if
      call default_bounds(n,lower,upper,lo,hi); call build_constraints(n,lo,hi,budget=.true.,c=c,cv=cv,d=d,dv=dv)
      mu=column_means(returns); r0=returns; if(dm) r0=returns-spread(mu,1,ns)
      if(present(min_return)) then
         if(present(mean_returns)) mu=mean_returns
         m=size(d,1); allocate(d2(m+1,n),dv2(m+1)); if(m>0) then; d2(1:m,:)=d; dv2(1:m)=dv; end if
         d2(m+1,:)=-mu; dv2(m+1)=-min_return; call move_alloc(d2,d); call move_alloc(dv2,dv)
      end if
      allocate(x0(n),grad(n),rw(ns)); x0=1.0_dp/real(n,dp); call project_feasible(x0,c,cv,d,dv,w,info); old=huge(1.0_dp)
      do it=1,mit
         rw=matmul(r0,w); grad=sum(r0*spread(sign(1.0_dp,rw),2,n),dim=1)/real(ns,dp); step=0.05_dp/sqrt(real(it,dp))
         x0=w-step*grad; call project_feasible(x0,c,cv,d,dv,w,info); obj=sum(abs(matmul(r0,w)))/real(ns,dp)
         if(abs(obj-old)<eps) exit; old=obj
      end do
      if(present(status)) status=merge(nmof_ok,nmof_numerical_failure,it<=mit)
   end subroutine minimum_mad

   subroutine default_bounds(n,lower,upper,lo,hi)
      integer,intent(in)::n
      real(dp),intent(in),optional::lower(:),upper(:)
      real(dp),allocatable,intent(out)::lo(:),hi(:)
      allocate(lo(n),hi(n)); lo=0.0_dp; hi=1.0_dp
      if(present(lower)) then; if(size(lower)==1) lo=lower(1); if(size(lower)==n) lo=lower; end if
      if(present(upper)) then; if(size(upper)==1) hi=upper(1); if(size(upper)==n) hi=upper; end if
   end subroutine default_bounds

   subroutine default_unbounded(n,lower,upper,lo,hi)
      integer,intent(in)::n
      real(dp),intent(in),optional::lower(:),upper(:)
      real(dp),allocatable,intent(out)::lo(:),hi(:)
      allocate(lo(n),hi(n)); lo=-huge(1.0_dp); hi=huge(1.0_dp)
      if(present(lower)) then; if(size(lower)==1) lo=lower(1); if(size(lower)==n) lo=lower; end if
      if(present(upper)) then; if(size(upper)==1) hi=upper(1); if(size(upper)==n) hi=upper; end if
   end subroutine default_unbounded

   subroutine build_constraints(n,lo,hi,groups,group_lower,group_upper,budget,c,cv,d,dv)
      integer,intent(in)::n
      real(dp),intent(in),optional::lo(:),hi(:),groups(:,:),group_lower(:),group_upper(:)
      logical,intent(in),optional::budget
      real(dp),allocatable,intent(out)::c(:,:),cv(:),d(:,:),dv(:)
      integer::m,i,k,ng
      logical::bud
      bud=.false.; if(present(budget)) bud=budget
      if(bud) then; allocate(c(1,n),cv(1)); c=1.0_dp; cv=1.0_dp; else; allocate(c(0,n),cv(0)); end if
      m=0
      if(present(hi)) m=m+count(hi<huge(1.0_dp)/4)
      if(present(lo)) m=m+count(lo>-huge(1.0_dp)/4)
      ng=0; if(present(groups)) ng=size(groups,1)
      if(present(group_upper)) m=m+size(group_upper)
      if(present(group_lower)) m=m+size(group_lower)
      allocate(d(m,n),dv(m)); k=0
      if(present(hi)) then; do i=1,n; if(hi(i)<huge(1.0_dp)/4) then; k=k+1; d(k,:)=0; d(k,i)=1; dv(k)=hi(i); end if; end do; end if
      if(present(lo)) then; do i=1,n; if(lo(i)>-huge(1.0_dp)/4) then; k=k+1; d(k,:)=0; d(k,i)=-1; dv(k)=-lo(i); end if; end do; end if
      if(present(group_upper).and.present(groups)) then; do i=1,size(group_upper); k=k+1; d(k,:)=groups(i,:); dv(k)=group_upper(i); end do; end if
      if(present(group_lower).and.present(groups)) then; do i=1,size(group_lower); k=k+1; d(k,:)=-groups(i,:); dv(k)=-group_lower(i); end do; end if
   end subroutine build_constraints

   subroutine regularize(q)
      real(dp),intent(inout)::q(:,:)
      integer::i
      do i=1,size(q,1); q(i,i)=q(i,i)+1e-10_dp*max(1.0_dp,maxval(abs(q))); end do
   end subroutine regularize

   pure function lowercase(s) result(t)
      character(len=*),intent(in)::s; character(len=len(s))::t; integer::i,c0
      do i=1,len(s); c0=iachar(s(i:i)); if(c0>=65.and.c0<=90) then; t(i:i)=achar(c0+32); else; t(i:i)=s(i:i); end if; end do
   end function lowercase
end module nmof_portfolio
